# MafA Discovery: Integrated Genomic Explorer
# ──────────────────────────────────────────
# Setup Instructions for External Users:
# 1. Ensure R is installed (https://cran.r-project.org/)
# 2. Open R or RStudio
# 3. Run the following lines to install dependencies:
#    install.packages(c("shiny", "data.table", "ggplot2", "plotly"))
#
# 4. REQUIRED FILES (Ensure these 5 files are in the same folder):
#    - app.R                             (This script)
#    - MafA_Peaks_with_SNPs_v3.csv       (Peak data)
#    - Master_DEG_Strain_Comparison_v3.csv (DEG data)
#    - mouse_genes_mm39_v3.csv           (Genomic backbone)
#    - B6_SJL_prioritized_protein_coding_SNPs.csv (SNP data)
#
# 5. TO LAUNCH:
#    Open this file in RStudio and click 'Run App', or run: shiny::runApp()

library(shiny)
library(data.table)
library(ggplot2)
library(plotly)

# ── 1. DATA ENGINE ────────────────────────────────────────────────────────
prepare_data <- function() {
  # Use local paths for portability
  chr_map <- data.table(
    Chr = paste0("Chr", c(1:19, "X", "Y")),
    Length = as.numeric(c(195154279, 181755017, 159745316, 156860686, 151754605, 149582041, 
                          144995196, 130127694, 124359700, 130530862, 121973369, 120092757, 
                          120885674, 125139656, 104073947, 98008968, 95294699, 90720763, 
                          61420004, 169476575, 91455967))
  )
  chr_map[, Offset := cumsum(data.table::shift(Length, fill = 0, type = "lag"))]
  chr_map[, Offset_Mbp := Offset / 1e6]
  chr_map[, Midpoint_Mbp := Offset_Mbp + (Length / 2e6)]
  
  # Load Required CSV Files
  master    <- fread("Master_DEG_Strain_Comparison_v3.csv")
  mafa      <- fread("MafA_Peaks_with_SNPs_v3.csv")
  all_genes <- fread("mouse_genes_mm39_v3.csv")
  snps      <- fread("B6_SJL_prioritized_protein_coding_SNPs.csv")
  
  # Format SNP data
  snps[, Chr := paste0("Chr", chr)]
  snps[, pos := as.numeric(pos)]
  snps[, phastCons_score := as.numeric(phastCons_score)]
  snps[is.na(phastCons_score), phastCons_score := 0]
  snps[, impact_simple := ifelse(any_high == TRUE | impact %like% "HIGH", "HIGH", "MODERATE")]
  
  # Add Global Position for Manhattan
  mafa[chr_map, on = "Chr", GlobalPos_Mbp := (Mid / 1e6) + i.Offset_Mbp]
  snps[chr_map, on = "Chr", GlobalPos_Mbp := (pos / 1e6) + i.Offset_Mbp]
  
  return(list(degs = master, mafa = mafa, genes = all_genes, chr_map = chr_map, snps = snps))
}

# ── 2. UI ──────────────────────────────────────────────────────────────────
ui <- fluidPage(
  tags$head(tags$style(HTML("
    .well-meta { background: #fffdf5; border: 1.5px solid gold; padding: 12px; margin-top: 15px; }
    .meta-title { font-weight: bold; font-size: 1.1em; border-bottom: 1px solid #ddd; margin-bottom: 8px; color: #856404; }
    .btn-rezoom { background-color: #007bff; color: white; font-weight: bold; width: 100%; margin-bottom: 10px; }
    .deg-item { margin-bottom: 4px; font-weight: bold; font-size: 0.92em; line-height: 1.2; }
    #hover_tooltip {
      position: absolute;
      pointer-events: none;
      background-color: rgba(255, 255, 255, 0.95);
      padding: 8px;
      border: 1px solid #ccc;
      border-radius: 4px;
      box-shadow: 2px 2px 5px rgba(0,0,0,0.2);
      font-size: 13px;
      z-index: 1000;
    }
  "))),
  
  titlePanel(span("MafA Discovery: Integrated Genomic Explorer", style="font-weight:bold;")),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      actionButton("reset_view", "↺ Reset to Genome-Wide", class = "btn-rezoom"),
      selectizeInput("search_gene", "Search Gene Symbol:", choices = NULL),
      selectInput("zoom_mode", "View Mode:", choices = c("Genome-Wide", "Chromosome", "Locus Zoom")),
      uiOutput("chr_selector_ui"),
      hr(),
      checkboxGroupInput("show_cat", "Visible Categories:", 
                         choices = c("Shared", "C57_Specific", "SJL_Specific", "Discordant"),
                         selected = c("Shared", "C57_Specific", "SJL_Specific", "Discordant")),
      numericInput("win_kb", "Locus Window (Kbp):", value = 500),
      sliderInput("min_score", "Min Peak Score:", 0, 10000, 0),
      hr(),
      checkboxInput("show_snps_main", "Show Coding SNPs on Main Plot", value = TRUE),
      checkboxGroupInput("snp_impact", "Coding SNP Impact:", 
                         choices = c("HIGH", "MODERATE"),
                         selected = c("HIGH", "MODERATE")),
      sliderInput("min_phastcons", "Min phastCons Score:", min = 0, max = 1, value = 0.7, step = 0.05),
      uiOutput("metadata_panel")
    ),
    mainPanel(
      width = 9,
      plotlyOutput("manhattan", height = "500px"),
      hr(),
      plotOutput("schematic", height = "350px")
    )
  )
)

# ── 3. SERVER ──────────────────────────────────────────────────────────────
server <- function(input, output, session) {
  d <- prepare_data()
  v <- reactiveValues(active_pk = NULL, active_snp = NULL, reset_trigger = 0, user_zoom = NULL)
  cat_colors <- c("Shared"="#228B22", "C57_Specific"="#0000CC", "SJL_Specific"="#CC0000", "Discordant"="#FF8C00", "Non-DE"="#D3D3D3")
  
  updateSelectizeInput(session, "search_gene", choices = sort(unique(d$genes$Symbol)), server = TRUE)
  output$chr_selector_ui <- renderUI({ selectInput("sel_chr", "Select Chromosome:", choices = d$chr_map$Chr) })
  
  observeEvent(input$reset_view, { 
    v$active_pk <- NULL
    v$active_snp <- NULL
    v$user_zoom <- NULL
    v$reset_trigger <- v$reset_trigger + 1 
    updateSelectInput(session, "zoom_mode", selected = "Genome-Wide") 
  })
  observeEvent(input$zoom_mode, { v$user_zoom <- NULL })
  observeEvent(input$sel_chr, { v$user_zoom <- NULL })

  # Search Implementation: Find nearest peak to searched gene
  observeEvent(input$search_gene, {
    req(input$search_gene != "")
    gene_row <- d$genes[Symbol == input$search_gene]
    if(nrow(gene_row) > 0) {
      gene_pos <- ifelse(gene_row$Strand[1] == 1, gene_row$Start[1], gene_row$End[1])
      # Find peaks on same Chr
      pks_chrom <- d$mafa[Chr == gene_row$Chr[1]]
      if(nrow(pks_chrom) > 0) {
        idx <- which.min(abs(pks_chrom$Mid - gene_pos))
        v$active_pk <- pks_chrom[idx]
        v$active_snp <- NULL
        updateSelectInput(session, "zoom_mode", selected = "Locus Zoom")
      }
    }
  })

  filtered_peaks <- reactive({
    req(input$show_cat)
    # 1. Filter by score
    pks <- d$mafa[`Peak Score` >= input$min_score]
    if(nrow(pks) == 0) return(NULL)
    
    # 2. Map to DEGs (Genomic Search via data.table nearest rolling join)
    deg_cols <- d$degs[, .(Chr, deg_start = as.numeric(Start), actual_deg_start = as.numeric(Start),
                           Category, Direction, Symbol, LFC_C57 = log2FC_C57, LFC_SJL = log2FC_SJL)]
    setkey(deg_cols, Chr, deg_start)
    
    pks_q <- copy(pks)
    pks_q[, query_pos := as.numeric(Mid)]
    res <- deg_cols[pks_q, on = .(Chr, deg_start = query_pos), roll = "nearest"]
    res <- res[abs(Mid - actual_deg_start) <= (input$win_kb * 1000) & !is.na(Category)]
    
    if(nrow(res) == 0) return(NULL)
    
    # Final filter by selected categories
    res[Category %in% input$show_cat]
  })
  
  filtered_snps <- reactive({
    req(input$snp_impact)
    min_pc <- if (is.null(input$min_phastcons)) 0 else input$min_phastcons
    d$snps[impact_simple %in% input$snp_impact & phastCons_score >= min_pc]
  })
  
  # Capture user's interactive zoom state from plotly
  observe({
    relayout <- event_data("plotly_relayout", source = "manhattan")
    if (!is.null(relayout)) {
      if (!is.null(relayout[["xaxis.range[0]"]])) {
        v$user_zoom <- c(relayout[["xaxis.range[0]"]], relayout[["xaxis.range[1]"]])
      }
      if (!is.null(relayout[["xaxis.autorange"]])) {
        v$user_zoom <- NULL
      }
    }
  })

  # Reactive Click Selection (Using Plotly events - for both peaks and coding SNPs)
  observe({
    event <- event_data("plotly_click", source = "manhattan")
    if(!is.null(event)) {
      df <- filtered_peaks()
      snps_df <- filtered_snps()
      
      # Check peak hits first
      hit_pk <- NULL
      if (!is.null(df) && nrow(df) > 0) {
        pk_cand <- df[abs(GlobalPos_Mbp - event$x) < 0.15]
        if (nrow(pk_cand) > 0) {
          idx <- which.min(abs(pk_cand$GlobalPos_Mbp - event$x))
          hit_pk <- pk_cand[idx]
        }
      }
      
      if (!is.null(hit_pk)) {
        v$active_pk <- hit_pk
        v$active_snp <- NULL
      } else if (!is.null(snps_df) && nrow(snps_df) > 0) {
        # Check SNP hits
        snp_cand <- snps_df[abs(GlobalPos_Mbp - event$x) < 0.15]
        if (nrow(snp_cand) > 0) {
          idx <- which.min(abs(snp_cand$GlobalPos_Mbp - event$x))
          clicked_snp <- snp_cand[idx]
          v$active_snp <- clicked_snp
          # Select nearest peak on the same chromosome
          pks_chr <- d$mafa[Chr == clicked_snp$Chr]
          if (nrow(pks_chr) > 0) {
            nearest_idx <- which.min(abs(pks_chr$Mid - clicked_snp$pos))
            v$active_pk <- pks_chr[nearest_idx]
          }
        }
      }
    }
  })

  output$manhattan <- renderPlotly({
    df <- filtered_peaks()
    if(is.null(df) || nrow(df) == 0) return(NULL)
    df <- as.data.frame(df)
    
    snps_df <- filtered_snps()
    
    # Auto-X Scale
    x_range <- if(input$zoom_mode == "Genome-Wide") {
      c(0, max(d$chr_map$Offset_Mbp + d$chr_map$Length/1e6))
    } else if(input$zoom_mode == "Chromosome") {
      req(input$sel_chr)
      if (!is.null(v$user_zoom)) {
        v$user_zoom
      } else {
        r <- d$chr_map[d$chr_map$Chr == input$sel_chr, ]; c(r$Offset_Mbp, r$Offset_Mbp + r$Length/1e6)
      }
    } else {
      req(v$active_pk); c(v$active_pk$GlobalPos_Mbp - 0.5, v$active_pk$GlobalPos_Mbp + 0.5)
    }
    
    # Auto-Zoom Y-Axis
    visible_pts <- df[df$GlobalPos_Mbp >= x_range[1] & df$GlobalPos_Mbp <= x_range[2], ]
    y_max <- if(nrow(visible_pts) > 0) max(visible_pts$`Peak Score`) * 1.15 else 500
    
    # Tooltip text for peaks
    df$text_label <- paste0("<b>Peak: ", df$Symbol, "</b>\nPeak Score: ", round(df$`Peak Score`, 0), "\nSNPs: ", df$variant_count)
    
    # Combined fill palette (DEGs + SNP Impact)
    all_cat_colors <- c(cat_colors, "HIGH" = "#CC00CC", "MODERATE" = "#DAA520")
    
    p <- ggplot() +
      geom_rect(data = d$chr_map, aes(xmin=Offset_Mbp, xmax=Offset_Mbp+Length/1e6, ymin=-Inf, ymax=Inf), 
                fill=rep(c("white", "#f9f9f9"), length.out=21), inherit.aes=FALSE) +
      geom_point(data=df, aes(x=GlobalPos_Mbp, y=`Peak Score`, fill=Category, shape=Direction, size=variant_count, text=text_label), 
                 color="black", stroke=0.3, alpha=0.6)
    
    show_snps <- isTRUE(input$show_snps_main)
    
    # Plot Coding SNPs scaled to secondary y-axis if enabled
    if (show_snps && !is.null(snps_df) && nrow(snps_df) > 0) {
      vis_snps <- snps_df[GlobalPos_Mbp >= x_range[1] & GlobalPos_Mbp <= x_range[2]]
      if (nrow(vis_snps) > 0) {
        vis_snps <- as.data.frame(vis_snps)
        vis_snps$y_scaled <- vis_snps$phastCons_score * y_max
        vis_snps$snp_size  <- 0.1 + 2.9 * vis_snps$phastCons_score
        vis_snps$snp_alpha <- 0.2 + 0.7 * vis_snps$phastCons_score
        vis_snps$text_label <- paste0(
          "<b>Coding SNP: </b>", vis_snps$gene_symbol_1,
          "\nConsequence: ", sub(";.*", "", vis_snps$csq),
          "\nAA Change: ", sub(";.*", "", vis_snps$aa_change),
          "\nphastCons: ", sprintf("%.2f", vis_snps$phastCons_score),
          "\nImpact: ", vis_snps$impact_simple
        )
        p <- p + geom_point(data=vis_snps, 
                            aes(x=GlobalPos_Mbp, y=y_scaled, fill=impact_simple, size=I(snp_size), alpha=I(snp_alpha), text=text_label), 
                            shape=23, color="black", stroke=0.2, inherit.aes=FALSE)
      }
    }
    
    p <- p +
      scale_fill_manual(name="Category / Impact", values=all_cat_colors) + 
      scale_shape_manual(name="Direction", values=c("UP"=24, "DOWN"=25)) +
      scale_size_continuous(range=c(1.9, 5.6), guide="none") + 
      theme_minimal(base_size = 14) + labs(x="Chromosome", y="MafA Peak Score") +
      theme(axis.line=element_line(color="black", linewidth=0.3), panel.grid=element_blank())
    
    y_scale_args <- list(limits=c(0, y_max), expand=expansion(mult=c(0.02, 0.05)))
    if (show_snps) {
      y_scale_args$sec.axis <- sec_axis(~ . / y_max, name = "phastCons Conservation Score (0-1)")
    }
    
    if (input$zoom_mode == "Genome-Wide") {
      p <- p + 
        scale_x_continuous(breaks=d$chr_map$Midpoint_Mbp, labels=gsub("Chr","",d$chr_map$Chr), expand=c(0,0)) +
        do.call(scale_y_continuous, y_scale_args)
    } else if (input$zoom_mode == "Chromosome") {
      chr_info <- d$chr_map[Chr == input$sel_chr]
      chr_len_mbp <- chr_info$Length / 1e6
      tick_interval <- if(chr_len_mbp > 120) 10 else 5
      all_breaks <- seq(0, floor(chr_len_mbp), by = tick_interval / 2)
      global_breaks <- all_breaks + chr_info$Offset_Mbp
      tick_labels <- ifelse(all_breaks %% tick_interval == 0, as.character(as.integer(all_breaks)), "")
      p <- p + 
        scale_x_continuous(limits=x_range, breaks=global_breaks, labels=tick_labels, expand=c(0,0)) +
        do.call(scale_y_continuous, y_scale_args) +
        labs(x = paste0(input$sel_chr, " (Mbp)"))
    } else {
      p <- p + 
        scale_x_continuous(limits=x_range, breaks=d$chr_map$Midpoint_Mbp, labels=gsub("Chr","",d$chr_map$Chr), expand=c(0,0)) +
        do.call(scale_y_continuous, y_scale_args)
    }
    
    if(!is.null(v$active_pk)) {
      p <- p + geom_point(data=v$active_pk, aes(x=GlobalPos_Mbp, y=`Peak Score`), shape=23, fill="gold", size=8, stroke=1.2, inherit.aes=FALSE)
    }
    if(show_snps && !is.null(v$active_snp)) {
      snp_hl_pos <- (v$active_snp$pos / 1e6) + d$chr_map[Chr == v$active_snp$Chr, Offset_Mbp]
      snp_hl_y   <- v$active_snp$phastCons_score * y_max
      snp_hl_df  <- data.frame(x = snp_hl_pos, y = snp_hl_y)
      p <- p + geom_point(data=snp_hl_df, aes(x=x, y=y), shape=23, fill="gold", size=8.5, stroke=1.2, inherit.aes=FALSE)
    }
    
    ggplotly(p, tooltip="text", source="manhattan") %>% 
      layout(uirevision = v$reset_trigger, 
             margin = list(t = 50),
             xaxis = list(ticks = "outside", ticklen = 5, tickcolor = "black"),
             yaxis = list(ticks = "outside", ticklen = 5, tickcolor = "black")) %>% 
      config(displayModeBar = FALSE)
  })

  output$schematic <- renderPlot({
    req(v$active_pk); pk <- v$active_pk; win <- input$win_kb * 1000
    
    window_genes <- d$genes[Chr == pk$Chr & End >= (pk$Mid - win) & Start <= (pk$Mid + win)]
    if(nrow(window_genes) == 0) return(NULL)
    
    plot_df <- merge(window_genes, d$degs[, .(GeneId, Category)], by="GeneId", all.x=TRUE)
    plot_df[is.na(Category), Category := "Non-DE"]
    plot_df[!Category %in% c(input$show_cat, "Non-DE"), Category := "Non-DE"]
    plot_df <- plot_df[Category != "Non-DE" | Type == "protein_coding"]
    if(nrow(plot_df) == 0) return(NULL)
    
    plot_df <- plot_df[order(Start)]
    plot_df[, y_lev := 1]
    if(nrow(plot_df) > 1) {
      y_rows   <- seq(1.2, 0, length.out = 4)
      last_end <- rep(-Inf, length(y_rows))
      buffer   <- win * 0.15
      for(i in seq_len(nrow(plot_df))) {
        placed <- FALSE
        for(r in seq_along(y_rows)) {
          if(plot_df$Start[i] > last_end[r] + buffer) {
            plot_df$y_lev[i] <- y_rows[r]
            last_end[r] <- plot_df$End[i]
            placed <- TRUE
            break
          }
        }
        if(!placed) plot_df$y_lev[i] <- min(y_rows)
      }
    }
    
    plot_df[, tss := ifelse(Strand == 1, Start, End)]
    plot_df[, tip := ifelse(Strand == 1, tss + (win * 0.025), tss - (win * 0.025))]
    rect_h <- ifelse(nrow(plot_df) > 15, 0.08, 0.12)
    
    # Map coding SNPs in the window
    snp_mapped <- data.table()
    if (!is.null(input$snp_impact) && length(input$snp_impact) > 0) {
      min_pc <- if (is.null(input$min_phastcons)) 0 else input$min_phastcons
      win_snps <- d$snps[Chr == pk$Chr & pos >= (pk$Mid - win) & pos <= (pk$Mid + win) &
                           impact_simple %in% input$snp_impact & phastCons_score >= min_pc]
      if (nrow(win_snps) > 0) {
        m1 <- merge(win_snps, plot_df[, .(GeneId, y_lev)], by.x = "gene_id_1", by.y = "GeneId", all.x = FALSE)
        m2 <- merge(win_snps[!gene_id_1 %in% plot_df$GeneId], plot_df[, .(Symbol, y_lev)], by.x = "gene_symbol_1", by.y = "Symbol", all.x = FALSE)
        snp_mapped <- rbind(m1, m2, fill = TRUE)
      }
    }
    
    p_schem <- ggplot(plot_df) +
      geom_rect(aes(xmin=pk$Start, xmax=pk$End, ymin=-Inf, ymax=Inf), fill="#FFF9C4", alpha=0.4) +
      geom_vline(xintercept = pk$Mid, color="gold4", linetype="dashed") +
      geom_rect(aes(xmin=Start, xmax=End, ymin=y_lev-rect_h, ymax=y_lev+rect_h, fill=Category), color="black", linewidth=0.25) +
      geom_segment(data=plot_df[Category != "Non-DE"], aes(x=tss, xend=tip, y=y_lev, yend=y_lev), 
                   arrow=arrow(length=unit(0.2, "cm"), type="closed"), color="black", linewidth=0.5)

    if (nrow(snp_mapped) > 0) {
      snp_low  <- snp_mapped[phastCons_score < 0.7]
      snp_high <- snp_mapped[phastCons_score >= 0.7]
      
      if (nrow(snp_low) > 0) {
        p_schem <- p_schem + 
          geom_point(data = snp_low, aes(x = pos, y = y_lev, color = impact_simple, size = phastCons_score),
                     shape = 5, stroke = 1.2, inherit.aes = FALSE)
      }
      if (nrow(snp_high) > 0) {
        p_schem <- p_schem + 
          geom_point(data = snp_high, aes(x = pos, y = y_lev, fill = impact_simple, size = phastCons_score),
                     shape = 23, color = "black", stroke = 0.5, inherit.aes = FALSE)
      }
      p_schem <- p_schem + scale_size_continuous(range = c(2.5, 6), limits = c(0, 1), guide = "none")
    }

    all_fill_colors  <- c(cat_colors, "HIGH" = "#CC00CC", "MODERATE" = "#DAA520")
    all_color_colors <- c(cat_colors, "HIGH" = "#CC00CC", "MODERATE" = "#DAA520")

    p_schem +
      geom_text(aes(x=pk$Mid, y=1.55, label=paste("MafA peak", sub(".*peak_", "", pk$PeakID))), 
                color="gold4", fontface="bold", size=6, inherit.aes=FALSE) +
      geom_text(aes(x = (Start + End)/2, y = y_lev + rect_h + 0.08, label = Symbol, color = Category), fontface = "bold", size = 6) +
      scale_fill_manual(values=all_fill_colors) + scale_color_manual(values=all_color_colors) +
      scale_x_continuous(labels = function(x) format(x/1e6, digits=5)) +
      labs(x=paste(pk$Chr, "(Mbp)"), y="") +
      theme_minimal() + 
      theme(
        axis.text.y=element_blank(), 
        panel.grid=element_blank(), 
        axis.line.x = element_line(),
        axis.text.x = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold")
      )
  })
  
  output$metadata_panel <- renderUI({
    req(v$active_pk); pk <- v$active_pk
    win_size <- input$win_kb * 1000
    local_degs <- d$degs[Chr == pk$Chr & abs(Start - pk$Mid) <= win_size]
    
    selected_impacts <- if (is.null(input$snp_impact)) c("HIGH", "MODERATE") else input$snp_impact
    min_pc <- if (is.null(input$min_phastcons)) 0 else input$min_phastcons
    
    local_snps <- d$snps[Chr == pk$Chr & abs(pos - pk$Mid) <= win_size &
                           impact_simple %in% selected_impacts & phastCons_score >= min_pc]
    
    wellPanel(class = "well-meta",
              div(class = "meta-title", paste("Peak", sub(".*peak_", "", pk$PeakID))),
              p(strong("Location: "), round(pk$Mid/1e6, 3), " Mbp"),
              p(strong("SNPs: "), pk$variant_count),
              hr(),
              p(strong("Locus DEGs & log2FC:")),
              if(nrow(local_degs) > 0) {
                tagList(lapply(1:nrow(local_degs), function(i) {
                  row <- local_degs[i]
                  fc_text <- if(row$Category %in% c("Shared", "Discordant")) {
                    paste0(" (C57:", round(row$log2FC_C57, 2), ", SJL:", round(row$log2FC_SJL, 2), ")")
                  } else if(row$Category == "C57_Specific") {
                    paste0(" (C57:", round(row$log2FC_C57, 2), ")")
                  } else {
                    paste0(" (SJL:", round(row$log2FC_SJL, 2), ")")
                  }
                  fc_text <- gsub("NA", "−", fc_text)
                  div(class="deg-item", style=paste0("color:", cat_colors[row$Category]), paste0("• ", row$Symbol, fc_text))
                }))
              } else { p("No DEGs in window.", style="font-style:italic") },
              hr(),
              p(strong("Coding SNPs in Locus:")),
              if (nrow(local_snps) > 0) {
                genes_with_snps <- unique(local_snps$gene_symbol_1[local_snps$gene_symbol_1 != ""])
                tagList(lapply(genes_with_snps, function(g_sym) {
                  g_snps <- local_snps[gene_symbol_1 == g_sym]
                  n_high <- sum(g_snps$impact_simple == "HIGH")
                  n_mod  <- sum(g_snps$impact_simple == "MODERATE")
                  max_pc <- max(g_snps$phastCons_score, na.rm = TRUE)
                  
                  counts_str <- c()
                  if (n_high > 0) counts_str <- c(counts_str, paste0(n_high, " HIGH"))
                  if (n_mod > 0)  counts_str <- c(counts_str, paste0(n_mod, " MODERATE"))
                  
                  gene_color <- if (n_high > 0) "#CC00CC" else "#DAA520"
                  
                  hi_snps <- g_snps[phastCons_score >= 0.7]
                  
                  div(style = "margin-bottom: 6px;",
                    div(class = "deg-item", style = paste0("color:", gene_color),
                        paste0("• ", g_sym, " (", paste(counts_str, collapse = ", "), " | max pCons: ", sprintf("%.2f", max_pc), ")")
                    ),
                    if (nrow(hi_snps) > 0) {
                      tagList(lapply(1:nrow(hi_snps), function(j) {
                        s <- hi_snps[j]
                        csq_clean <- sub(";.*", "", s$csq)
                        aa_clean  <- sub(";.*", "", s$aa_change)
                        div(style = "font-size: 0.82em; font-family: monospace; color: #444; margin-left: 10px;",
                            paste0("★ pos:", s$pos, " | ", csq_clean, " | aa:", aa_clean, " | pCons:", sprintf("%.2f", s$phastCons_score))
                        )
                      }))
                    }
                  )
                }))
              } else { p("No coding SNPs in window.", style="font-style:italic") }
    )
  })
}

shinyApp(ui, server)
