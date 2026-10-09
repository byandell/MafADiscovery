# MafA Discovery: Integrated Genomic Explorer
# ──────────────────────────────────────────
# Setup Instructions for External Users:
# 1. Ensure R is installed (https://cran.r-project.org/)
# 2. Open R or RStudio
# 3. Run the following lines to install dependencies:
#    install.packages(c("shiny", "data.table", "ggplot2", "plotly"))
#
# 4. REQUIRED FILES (Ensure these 8 files are in the same folder):
#    - app.R                                           (This script)
#    - MafA_Peaks_with_SNPs_v3.csv                     (Peak data)
#    - Master_DEG_Strain_Comparison_v3.csv             (DEG data)
#    - mouse_genes_mm39_v3.csv                         (Genomic backbone)
#    - B6_SJL_prioritized_protein_coding_SNPs.csv      (SNP data)
#    - Top_glycemic_QTL_for_sex_additive_analysis.csv  (F2 Glycemic QTL data)
#    - interpretation_guide.md                         (Interpretation modal guide)
#    - developer_guide.md                              (Developer architecture modal)
#
# 5. TO LAUNCH:
#    Open this file in RStudio and click 'Run App', or run: shiny::runApp()

library(shiny)
library(data.table)
library(ggplot2)
library(plotly)

# Serve docs directory for local documentation links if available
docs_path <- if (dir.exists("docs")) "docs" else if (dir.exists("../docs")) "../docs" else NULL
if (!is.null(docs_path)) {
  shiny::addResourcePath("docs", normalizePath(docs_path))
}

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
  qtls      <- fread("Top_glycemic_QTL_for_sex_additive_analysis.csv")
  
  # Format SNP data
  snps[, Chr := paste0("Chr", chr)]
  snps[, pos := as.numeric(pos)]
  snps[, phastCons_score := as.numeric(phastCons_score)]
  snps[is.na(phastCons_score), phastCons_score := 0]
  snps[, impact_simple := ifelse(any_high == TRUE | impact %like% "HIGH", "HIGH", "MODERATE")]
  
  # Format QTL data
  qtls[, Chr := paste0("Chr", Chr)]
  qtls[, pos := as.numeric(pos)]
  qtls[, ci.low := as.numeric(ci.low)]
  qtls[, ci.high := as.numeric(ci.high)]
  qtls[, lod := as.numeric(lod)]
  qtls[, qtl_id := paste0(trait, " @ ", Chr, ":", round(pos, 1), " Mb (LOD ", round(lod, 1), ")")]
  
  # Add Global Position for Manhattan
  mafa[chr_map, on = "Chr", GlobalPos_Mbp := (Mid / 1e6) + i.Offset_Mbp]
  snps[chr_map, on = "Chr", GlobalPos_Mbp := (pos / 1e6) + i.Offset_Mbp]
  qtls[chr_map, on = "Chr", `:=`(
    GlobalPos_Mbp     = pos + i.Offset_Mbp,
    GlobalCI_Low_Mbp  = ci.low + i.Offset_Mbp,
    GlobalCI_High_Mbp = ci.high + i.Offset_Mbp
  )]
  
  return(list(degs = master, mafa = mafa, genes = all_genes, chr_map = chr_map, snps = snps, qtls = qtls))
}

# ── 2. UI ──────────────────────────────────────────────────────────────────
# ── 2. UI ──────────────────────────────────────────────────────────────────
ui <- fluidPage(
  title = "MafA Discovery: Integrated Genomic Explorer",
  tags$head(
    tags$title("MafA Discovery: Integrated Genomic Explorer"),
    tags$link(rel = "icon", href = "data:image/svg+xml,<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 100 100'><text y='.9em' font-size='90'>🧬</text></svg>"),
    tags$style(HTML("
      body { background-color: #f8f9fa; font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif; }
      .well-nav { background: #ffffff; border: 1px solid #d0d7de; border-radius: 8px; padding: 12px 18px; margin-top: 10px; margin-bottom: 16px; box-shadow: 0 1px 3px rgba(0,0,0,0.05); }
      .well-ctrl { background: #ffffff; border: 1px solid #d0d7de; border-radius: 8px; padding: 14px 16px; margin-bottom: 16px; box-shadow: 0 1px 3px rgba(0,0,0,0.04); font-size: 0.92em; }
      .ctrl-header { display: flex; justify-content: space-between; align-items: center; border-bottom: 1.5px solid #edf0f2; padding-bottom: 8px; margin-bottom: 10px; }
      .ctrl-header .checkbox { margin: 0; }
      .panel-container { background: #ffffff; border: 1px solid #d0d7de; border-radius: 8px; padding: 10px; margin-bottom: 16px; box-shadow: 0 1px 3px rgba(0,0,0,0.04); }
      .well-meta { background: #fffdf5; border: 1.5px solid gold; padding: 16px; margin-bottom: 20px; border-radius: 8px; box-shadow: 0 1px 3px rgba(0,0,0,0.04); }
      .meta-title { font-weight: bold; font-size: 1.18em; color: #856404; }
      .btn-rezoom { background-color: #007bff; color: white; font-weight: bold; }
      .btn-rezoom:hover { background-color: #0056b3; color: white; }
      .btn-help { background-color: #17a2b8; color: white; font-weight: bold; }
      .btn-help:hover { background-color: #117a8b; color: white; }
      .btn-qtl { background-color: #6f42c1; color: white; font-weight: bold; }
      .btn-qtl:hover { background-color: #59359a; color: white; }
      .btn-devguide { background-color: #495057; color: white; font-weight: bold; }
      .btn-devguide:hover { background-color: #343a40; color: white; }
      .deg-item { margin-bottom: 4px; font-weight: bold; font-size: 0.92em; line-height: 1.25; }
      .modal-markdown { line-height: 1.6; font-size: 0.96em; }
      .modal-markdown h1 { font-size: 1.45em; font-weight: bold; margin-bottom: 14px; color: #222; }
      .modal-markdown h2 { font-size: 1.25em; font-weight: bold; margin-top: 18px; color: #333; }
      .modal-markdown h3 { font-size: 1.1em; font-weight: bold; margin-top: 14px; color: #444; }
      .modal-markdown a { color: #007bff; text-decoration: underline; font-weight: bold; }
      .modal-markdown a:hover { color: #0056b3; }
      .modal-markdown ul, .modal-markdown ol { padding-left: 24px; margin-bottom: 12px; }
      .modal-markdown li { margin-bottom: 6px; }
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
    ")),
    tags$script(HTML("
      document.title = 'MafA Discovery: Integrated Genomic Explorer';
      window.name = 'mafa_app';
      $(document).on('shown.bs.modal', function () {
        $('.modal-markdown a').attr('target', '_blank').attr('rel', 'noopener noreferrer');
      });
      function resolveDocsUrl() {
        if (window.location.pathname.indexOf('/app_') !== -1 || window.location.hostname.includes('github.io')) {
          var base = window.location.pathname.replace(/\\/app_[^\\/]+.*$/, '').replace(/\\/index\\.html$/, '').replace(/\\/$/, '');
          return window.location.origin + base + '/docs';
        }
        return 'docs';
      }
      $(document).on('click', '#btn-devguide-link', function(e) {
        $(this).attr('href', resolveDocsUrl());
      });
    "))
  ),
  
  # ── Top Global Navigation Bar ──────────────────────────────────────────
  div(class = "well-nav",
    fluidRow(
      column(width = 2,
        div(style = "padding-top: 4px;",
          span("🧬 MafA Discovery", style = "font-weight: 800; font-size: 1.28em; color: #1976d2; display: block; line-height: 1.15;"),
          span("Genomic Explorer", style = "font-size: 0.84em; color: #666; font-weight: 500;")
        )
      ),
      column(width = 3,
        selectizeInput("search_gene", "Search Gene Symbol:", choices = NULL, width = "100%")
      ),
      column(width = 2,
        selectInput("zoom_mode", "View Mode:", choices = c("Genome-Wide", "Chromosome", "QTL Region"), width = "100%")
      ),
      column(width = 2,
        uiOutput("context_selector_ui")
      ),
      column(width = 3,
        div(style = "display: flex; gap: 6px; flex-wrap: wrap; justify-content: flex-end; margin-top: 24px;",
          actionButton("show_qtl_table", "📊 QTLs", class = "btn btn-sm btn-qtl"),
          actionButton("show_help", "ℹ️ Guide", class = "btn btn-sm btn-help"),
          tags$a(id = "btn-devguide-link", href = "docs", target = "_blank", rel = "opener",
                 onclick = "this.href = resolveDocsUrl();",
                 class = "btn btn-sm btn-devguide", style = "text-decoration: none; display: inline-flex; align-items: center;", "📖 Dev Guide ↗"),
          actionButton("reset_view", "↺ Reset", class = "btn btn-sm btn-rezoom")
        )
      )
    )
  ),
  
  # ── Macro Manhattan Plot Row ───────────────────────────────────────────
  fluidRow(
    column(width = 9,
      div(class = "panel-container",
        plotlyOutput("manhattan", height = "480px")
      )
    ),
    column(width = 3,
      wellPanel(class = "well-ctrl",
        div(class = "ctrl-header",
          span(strong("Manhattan Controls"), style = "color: #1976d2; font-size: 1.05em;"),
          div(style = "margin: 0;",
            checkboxInput("hide_manhattan_legend", "Hide Legend", value = FALSE)
          )
        ),
        checkboxGroupInput("show_cat", "Visible Peak Categories:", 
                           choices = c("Shared", "C57_Specific", "SJL_Specific", "Discordant"),
                           selected = c("C57_Specific", "SJL_Specific", "Discordant")),
        selectInput("min_score", "Min Peak Score:", 
                    choices = c("All (0)" = 0, "100" = 100, "200" = 200, 
                                "500" = 500, "1,000" = 1000, "2,000" = 2000, 
                                "5,000" = 5000), 
                    selected = 0),
        hr(style = "margin: 8px 0;"),
        checkboxInput("show_snps_main", "Show Coding SNPs on Main Plot", value = TRUE),
        conditionalPanel(
          condition = "input.show_snps_main == true",
          checkboxGroupInput("snp_impact", "Coding SNP Impact:", 
                             choices = c("HIGH", "MODERATE"),
                             selected = c("HIGH", "MODERATE")),
          sliderInput("min_phastcons", "Min phastCons Score:", min = 0, max = 1, value = 0.7, step = 0.05)
        )
      )
    )
  ),
  
  # ── Micro Locus Schematic Row ──────────────────────────────────────────
  fluidRow(
    column(width = 9,
      div(class = "panel-container",
        plotOutput("schematic", height = "340px")
      )
    ),
    column(width = 3,
      wellPanel(class = "well-ctrl",
        div(class = "ctrl-header",
          span(strong("Locus Schematic Controls"), style = "color: #856404; font-size: 1.05em;"),
          div(style = "margin: 0;",
            checkboxInput("hide_schematic_legend", "Hide Legend", value = FALSE)
          )
        ),
        selectInput("win_kb", "Locus Window Size:", 
                    choices = c("20 kb" = 20, "50 kb" = 50, "100 kb" = 100, 
                                "200 kb" = 200, "500 kb" = 500, "1,000 kb (1 Mb)" = 1000, 
                                "2,000 kb (2 Mb)" = 2000), 
                    selected = 500),
        div(style = "font-size: 0.85em; color: #555; margin-top: 12px; line-height: 1.45; background: #fffdf5; padding: 10px; border-radius: 4px; border: 1px solid #fae8a4;",
          p(style = "margin-bottom: 5px;", strong("Micro-Architecture Guide:")),
          p(style = "margin-bottom: 4px;", "• ", strong("TSS arrows:"), " Transcription start site & direction."),
          p(style = "margin-bottom: 4px;", "• ", strong("Gold line:"), " Active peak center."),
          p(style = "margin-bottom: 4px;", "• ", strong("Diamonds:"), " High/moderate coding SNPs."),
          p(style = "margin-bottom: 0;", "• Click any point in the Manhattan plot to reposition.")
        )
      )
    )
  ),
  
  # ── Peak Metadata Row ──────────────────────────────────────────────────
  fluidRow(
    column(width = 12,
      uiOutput("metadata_panel")
    )
  )
)

# ── 3. SERVER ──────────────────────────────────────────────────────────────
server <- function(input, output, session) {
  d <- prepare_data()
  v <- reactiveValues(
    active_pk   = NULL, 
    active_snp  = NULL, 
    active_qtl  = d$qtls[1], 
    last_gene   = NULL, 
    current_chr = "Chr1", 
    reset_trigger = 0, 
    user_zoom   = NULL
  )
  cat_colors <- c("Shared"="#228B22", "C57_Specific"="#0000CC", "SJL_Specific"="#CC0000", "Discordant"="#FF8C00", "Non-DE"="#D3D3D3")
  
  updateSelectizeInput(session, "search_gene", choices = c("", sort(unique(d$genes$Symbol))), server = TRUE)
  
  # Discrete input accessors with fallbacks
  get_win_kb <- reactive({
    if (is.null(input$win_kb)) 500 else as.numeric(input$win_kb)
  })
  
  get_min_score <- reactive({
    if (is.null(input$min_score)) 0 else as.numeric(input$min_score)
  })
  
  # Dynamic context selector for top nav bar (Chromosome, QTL, or indicator)
  output$context_selector_ui <- renderUI({
    if (input$zoom_mode == "Chromosome") {
      selectInput("sel_chr", "Select Chromosome:", choices = d$chr_map$Chr, selected = v$current_chr, width = "100%")
    } else if (input$zoom_mode == "QTL Region") {
      sel_val <- if (!is.null(v$active_qtl)) v$active_qtl$qtl_id else d$qtls$qtl_id[1]
      selectInput("sel_qtl", "Select F2 Glycemic QTL:", choices = d$qtls$qtl_id, selected = sel_val, width = "100%")
    } else if (input$zoom_mode == "Locus Zoom") {
      div(style = "padding-top: 25px;",
        span(class = "badge", style = "background-color: #856404; font-size: 0.85em; padding: 6px 10px;",
             if (!is.null(v$active_pk)) paste("Locus:", v$active_pk$Symbol) else "Locus Zoom")
      )
    } else {
      div(style = "padding-top: 25px; color: #6c757d; font-size: 0.88em; font-style: italic;",
          "Genome-wide view")
    }
  })
  
  # Centralized state reset function
  reset_to_default_state <- function() {
    v$active_pk   <- NULL
    v$active_snp  <- NULL
    v$last_gene   <- NULL
    v$current_chr <- "Chr1"
    v$active_qtl  <- d$qtls[1]
    v$user_zoom   <- NULL
    v$reset_trigger <- v$reset_trigger + 1 
    
    # Completely clear gene search selectize
    updateSelectizeInput(session, "search_gene", choices = c("", sort(unique(d$genes$Symbol))), selected = "", server = TRUE)
    
    # Reset View Mode choices and selection to Genome-Wide
    updateSelectInput(session, "zoom_mode", 
                      choices = c("Genome-Wide", "Chromosome", "QTL Region"), 
                      selected = "Genome-Wide")
    
    # Reset Category Checkboxes (Shared deselected)
    updateCheckboxGroupInput(session, "show_cat", 
                             selected = c("C57_Specific", "SJL_Specific", "Discordant"))
    
    # Reset Min Peak Score
    updateSelectInput(session, "min_score", selected = "0")
    
    # Reset Coding SNP controls
    updateCheckboxInput(session, "show_snps_main", value = TRUE)
    updateCheckboxGroupInput(session, "snp_impact", selected = c("HIGH", "MODERATE"))
    updateSliderInput(session, "min_phastcons", value = 0.7)
    
    # Reset Locus Window size
    if (!is.null(input$win_kb) && input$win_kb != "500") {
      updateSelectInput(session, "win_kb", selected = "500")
    }
    
    # Reset legend toggles
    if (isTRUE(input$hide_manhattan_legend)) {
      updateCheckboxInput(session, "hide_manhattan_legend", value = FALSE)
    }
    if (isTRUE(input$hide_schematic_legend)) {
      updateCheckboxInput(session, "hide_schematic_legend", value = FALSE)
    }
  }
  
  # Reset view action button
  observeEvent(input$reset_view, { 
    reset_to_default_state()
  })
  
  # Auto-reset on "Genome-Wide" view mode selection
  observeEvent(input$zoom_mode, { 
    v$user_zoom <- NULL 
    if (input$zoom_mode == "Genome-Wide") {
      has_custom_state <- !is.null(v$active_pk) || !is.null(v$last_gene) || 
        (!is.null(input$search_gene) && nzchar(input$search_gene)) ||
        !is.null(v$active_snp) || 
        !identical(sort(input$show_cat), sort(c("C57_Specific", "SJL_Specific", "Discordant"))) ||
        (!is.null(input$min_score) && input$min_score != "0") || 
        (!is.null(input$show_snps_main) && !isTRUE(input$show_snps_main)) || 
        (!is.null(input$min_phastcons) && input$min_phastcons != 0.7) ||
        isTRUE(input$hide_manhattan_legend) || isTRUE(input$hide_schematic_legend)
      
      if (has_custom_state) {
        reset_to_default_state()
      }
    } else if (input$zoom_mode == "QTL Region" && !is.null(v$active_qtl)) {
      v$last_gene <- NULL
      if (!is.null(input$search_gene) && nzchar(input$search_gene)) {
        updateSelectizeInput(session, "search_gene", choices = c("", sort(unique(d$genes$Symbol))), selected = "", server = TRUE)
      }
      v$current_chr <- v$active_qtl$Chr
      local_pks <- d$mafa[Chr == v$active_qtl$Chr & Mid >= v$active_qtl$ci.low * 1e6 & Mid <= v$active_qtl$ci.high * 1e6]
      if (nrow(local_pks) > 0) {
        v$active_pk <- local_pks[which.max(`Peak Score`)]
      }
    }
  })
  
  observeEvent(input$sel_chr, { 
    v$user_zoom <- NULL
    if (!is.null(input$sel_chr) && nzchar(input$sel_chr)) {
      v$current_chr <- input$sel_chr
    }
  })

  # QTL Selection observer
  observeEvent(input$sel_qtl, {
    req(input$sel_qtl)
    q_row <- d$qtls[qtl_id == input$sel_qtl]
    if (nrow(q_row) > 0) {
      v$active_qtl   <- q_row[1]
      v$current_chr  <- q_row$Chr[1]
      v$user_zoom    <- NULL
      v$last_gene    <- NULL
      if (!is.null(input$search_gene) && nzchar(input$search_gene)) {
        updateSelectizeInput(session, "search_gene", choices = c("", sort(unique(d$genes$Symbol))), selected = "", server = TRUE)
      }
      local_pks <- d$mafa[Chr == q_row$Chr[1] & Mid >= q_row$ci.low * 1e6 & Mid <= q_row$ci.high * 1e6]
      if (nrow(local_pks) > 0) {
        v$active_pk <- local_pks[which.max(`Peak Score`)]
      }
    }
  })

  # Search Implementation: Find nearest peak to searched gene, or revert when deselected
  observeEvent(input$search_gene, {
    if (!is.null(input$search_gene) && nzchar(input$search_gene)) {
      gene_row <- d$genes[Symbol == input$search_gene]
      if(nrow(gene_row) > 0) {
        gene_pos <- ifelse(gene_row$Strand[1] == 1, gene_row$Start[1], gene_row$End[1])
        # Find peaks on same Chr
        pks_chrom <- d$mafa[Chr == gene_row$Chr[1]]
        if(nrow(pks_chrom) > 0) {
          idx <- which.min(abs(pks_chrom$Mid - gene_pos))
          v$active_pk   <- pks_chrom[idx]
          v$active_snp  <- NULL
          v$last_gene   <- input$search_gene
          v$current_chr <- gene_row$Chr[1]
          updateSelectInput(session, "zoom_mode", 
                            choices = c("Genome-Wide", "Chromosome", "QTL Region", "Locus Zoom"), 
                            selected = "Locus Zoom")
        }
      }
    } else {
      # Deselected / cleared: revert to whole genome only if in Locus Zoom mode from search
      if (identical(input$zoom_mode, "Locus Zoom") && !is.null(v$last_gene)) {
        reset_to_default_state()
      } else {
        v$last_gene <- NULL
      }
    }
  }, ignoreInit = TRUE)

  # Interactive QTL Reference Table Modal
  observeEvent(input$show_qtl_table, {
    showModal(modalDialog(
      title = span(strong("Top Glycemic QTLs (F2 Study - Sex Additive Analysis)"), style = "color: #6f42c1;"),
      div(
        p("Select a QTL below to jump directly to its confidence interval boundaries on the Manhattan plot:"),
        div(style = "display: flex; gap: 10px; align-items: flex-end; margin-bottom: 15px;",
          div(style = "flex: 1;",
            selectInput("modal_qtl_select", "Select QTL:", choices = d$qtls$qtl_id, 
                        selected = if(!is.null(v$active_qtl)) v$active_qtl$qtl_id else d$qtls$qtl_id[1],
                        width = "100%")
          ),
          actionButton("jump_qtl_btn", "Zoom to QTL", class = "btn-qtl", 
                       style = "margin-bottom: 15px; height: 38px;")
        ),
        hr(),
        div(style = "max-height: 400px; overflow-y: auto;",
          tableOutput("qtl_summary_table")
        )
      ),
      size = "l",
      easyClose = TRUE,
      footer = modalButton("Close")
    ))
  })

  observeEvent(input$jump_qtl_btn, {
    req(input$modal_qtl_select)
    q_row <- d$qtls[qtl_id == input$modal_qtl_select]
    if (nrow(q_row) > 0) {
      v$active_qtl   <- q_row[1]
      v$current_chr  <- q_row$Chr[1]
      v$user_zoom    <- NULL
      v$last_gene    <- NULL
      if (!is.null(input$search_gene) && nzchar(input$search_gene)) {
        updateSelectizeInput(session, "search_gene", choices = c("", sort(unique(d$genes$Symbol))), selected = "", server = TRUE)
      }
      local_pks <- d$mafa[Chr == q_row$Chr[1] & Mid >= q_row$ci.low * 1e6 & Mid <= q_row$ci.high * 1e6]
      if (nrow(local_pks) > 0) {
        v$active_pk <- local_pks[which.max(`Peak Score`)]
      }
      updateSelectInput(session, "zoom_mode", 
                        choices = c("Genome-Wide", "Chromosome", "QTL Region", if (!is.null(v$active_pk)) "Locus Zoom"), 
                        selected = "QTL Region")
      removeModal()
    }
  })

  output$qtl_summary_table <- renderTable({
    tab <- copy(d$qtls)
    tab[, `:=`(
      Trait = trait,
      Marker = marker,
      Chr = Chr,
      `Peak (Mb)` = sprintf("%.2f", pos),
      `95% CI (Mb)` = paste0(sprintf("%.2f", ci.low), " – ", sprintf("%.2f", ci.high)),
      LOD = sprintf("%.2f", lod),
      `BB Effect` = sprintf("%.3f", BB_effect),
      `BS Effect` = sprintf("%.3f", BS_effect),
      `SS Effect` = sprintf("%.3f", SS_effect)
    )]
    tab[, .(Trait, Marker, Chr, `Peak (Mb)`, `95% CI (Mb)`, LOD, `BB Effect`, `BS Effect`, `SS Effect`)]
  }, striped = TRUE, hover = TRUE, bordered = TRUE, spacing = "s")

  # Helper to load and render markdown files inside modals
  render_markdown_file <- function(filename) {
    filepath <- if (file.exists(filename)) {
      filename
    } else if (file.exists(file.path("SHINY_APP", filename))) {
      file.path("SHINY_APP", filename)
    } else {
      NULL
    }
    
    if (!is.null(filepath)) {
      content <- paste(readLines(filepath, encoding = "UTF-8", warn = FALSE), collapse = "\n")
      div(class = "modal-markdown",
        if (exists("markdown", where = asNamespace("shiny"), mode = "function")) {
          shiny::markdown(content)
        } else if (requireNamespace("markdown", quietly = TRUE)) {
          HTML(markdown::markdownToHTML(text = content, fragment.only = TRUE))
        } else {
          tags$pre(content)
        }
      )
    } else {
      tags$p(paste("Documentation file not found:", filename), style = "color: red;")
    }
  }

  # Interpretation Modal Guide
  observeEvent(input$show_help, {
    showModal(modalDialog(
      title = span(strong("MafA Discovery: How to Interpret Panels & Navigate"), style = "color: #17a2b8;"),
      render_markdown_file("interpretation_guide.md"),
      size = "l",
      easyClose = TRUE,
      footer = modalButton("Close")
    ))
  })

  # Developer Guide Modal
  observeEvent(input$show_dev_guide, {
    showModal(modalDialog(
      title = span(strong("MafA Discovery: Developer Guide & Architecture"), style = "color: #343a40;"),
      render_markdown_file("developer_guide.md"),
      size = "l",
      easyClose = TRUE,
      footer = modalButton("Close")
    ))
  })

  filtered_peaks <- reactive({
    req(input$show_cat)
    # 1. Filter by score
    pks <- d$mafa[`Peak Score` >= get_min_score()]
    if(nrow(pks) == 0) return(NULL)
    
    # 2. Map to DEGs (Genomic Search via data.table nearest rolling join)
    deg_cols <- d$degs[, .(Chr, deg_start = as.numeric(Start), actual_deg_start = as.numeric(Start),
                           Category, Direction, Symbol, LFC_C57 = log2FC_C57, LFC_SJL = log2FC_SJL)]
    setkey(deg_cols, Chr, deg_start)
    
    pks_q <- copy(pks)
    pks_q[, query_pos := as.numeric(Mid)]
    res <- deg_cols[pks_q, on = .(Chr, deg_start = query_pos), roll = "nearest"]
    res <- res[abs(Mid - actual_deg_start) <= (get_win_kb() * 1000) & !is.na(Category)]
    
    if(nrow(res) == 0) return(NULL)
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
        v$active_pk   <- hit_pk
        v$active_snp  <- NULL
        v$current_chr <- hit_pk$Chr
        updateSelectInput(session, "zoom_mode", 
                          choices = c("Genome-Wide", "Chromosome", "QTL Region", "Locus Zoom"),
                          selected = input$zoom_mode)
      } else if (!is.null(snps_df) && nrow(snps_df) > 0) {
        # Check SNP hits
        snp_cand <- snps_df[abs(GlobalPos_Mbp - event$x) < 0.15]
        if (nrow(snp_cand) > 0) {
          idx <- which.min(abs(snp_cand$GlobalPos_Mbp - event$x))
          clicked_snp <- snp_cand[idx]
          v$active_snp  <- clicked_snp
          v$current_chr <- clicked_snp$Chr
          # Select nearest peak on the same chromosome
          pks_chr <- d$mafa[Chr == clicked_snp$Chr]
          if (nrow(pks_chr) > 0) {
            nearest_idx <- which.min(abs(pks_chr$Mid - clicked_snp$pos))
            v$active_pk <- pks_chr[nearest_idx]
          }
          updateSelectInput(session, "zoom_mode", 
                            choices = c("Genome-Wide", "Chromosome", "QTL Region", "Locus Zoom"),
                            selected = input$zoom_mode)
        }
      }
    }
  })

  output$manhattan <- renderPlotly({
    df <- filtered_peaks()
    if(is.null(df) || nrow(df) == 0) return(NULL)
    df <- as.data.frame(df)
    
    snps_df <- filtered_snps()
    
    cur_chr <- if(!is.null(input$sel_chr) && nzchar(input$sel_chr)) input$sel_chr else v$current_chr
    
    # Auto-X Scale
    x_range <- if(input$zoom_mode == "Genome-Wide") {
      c(0, max(d$chr_map$Offset_Mbp + d$chr_map$Length/1e6))
    } else if(input$zoom_mode == "Chromosome") {
      if (!is.null(v$user_zoom)) {
        v$user_zoom
      } else {
        r <- d$chr_map[d$chr_map$Chr == cur_chr, ]
        c(r$Offset_Mbp, r$Offset_Mbp + r$Length/1e6)
      }
    } else if(input$zoom_mode == "QTL Region") {
      req(v$active_qtl)
      q <- v$active_qtl
      chr_offset <- d$chr_map[Chr == q$Chr, Offset_Mbp]
      ci_span <- q$ci.high - q$ci.low
      pad <- max(0.5, ci_span * 0.03)
      if (!is.null(v$user_zoom)) {
        v$user_zoom
      } else {
        c(chr_offset + q$ci.low - pad, chr_offset + q$ci.high + pad)
      }
    } else {
      req(v$active_pk)
      win_mbp <- get_win_kb() / 1000
      if (!is.null(v$user_zoom)) {
        v$user_zoom
      } else {
        c(v$active_pk$GlobalPos_Mbp - win_mbp, v$active_pk$GlobalPos_Mbp + win_mbp)
      }
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
                fill=rep(c("white", "#f9f9f9"), length.out=21), inherit.aes=FALSE)
    
    # Highlight active QTL interval if in QTL Region or matching Chromosome
    if (!is.null(v$active_qtl)) {
      q <- v$active_qtl
      if (input$zoom_mode == "QTL Region" || (input$zoom_mode == "Chromosome" && cur_chr == q$Chr)) {
        p <- p + 
          annotate("rect", xmin = q$GlobalCI_Low_Mbp, xmax = q$GlobalCI_High_Mbp, ymin = -Inf, ymax = Inf,
                   fill = "#007bff", alpha = 0.08) +
          geom_vline(xintercept = q$GlobalPos_Mbp, color = "#007bff", linetype = "dashed", linewidth = 0.6)
      }
    }
    
    p <- p + geom_point(data=df, aes(x=GlobalPos_Mbp, y=`Peak Score`, fill=Category, shape=Direction, size=variant_count, text=text_label), 
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
      chr_info <- d$chr_map[Chr == cur_chr]
      chr_len_mbp <- chr_info$Length / 1e6
      tick_interval <- if(chr_len_mbp > 120) 10 else 5
      all_breaks <- seq(0, floor(chr_len_mbp), by = tick_interval / 2)
      global_breaks <- all_breaks + chr_info$Offset_Mbp
      tick_labels <- ifelse(all_breaks %% tick_interval == 0, as.character(as.integer(all_breaks)), "")
      p <- p + 
        scale_x_continuous(limits=x_range, breaks=global_breaks, labels=tick_labels, expand=c(0,0)) +
        do.call(scale_y_continuous, y_scale_args) +
        labs(x = paste0(cur_chr, " (Mbp)"))
    } else if (input$zoom_mode == "QTL Region") {
      req(v$active_qtl)
      q <- v$active_qtl
      chr_offset <- d$chr_map[Chr == q$Chr, Offset_Mbp]
      local_min <- max(0, x_range[1] - chr_offset)
      local_max <- min(d$chr_map[Chr == q$Chr, Length / 1e6], x_range[2] - chr_offset)
      local_breaks <- pretty(c(local_min, local_max), n = 7)
      local_breaks <- local_breaks[local_breaks >= local_min & local_breaks <= local_max]
      global_breaks <- local_breaks + chr_offset
      diff_mbp <- local_max - local_min
      dec <- if (diff_mbp <= 1) 2 else 1
      local_labels <- sprintf(paste0("%.", dec, "f"), local_breaks)
      p <- p + 
        scale_x_continuous(limits=x_range, breaks=global_breaks, labels=local_labels, expand=c(0,0)) +
        do.call(scale_y_continuous, y_scale_args) +
        labs(x = paste0(q$Chr, " (Mbp) — ", q$trait, " QTL [95% CI: ", round(q$ci.low, 1), "–", round(q$ci.high, 1), " Mb]"))
    } else {
      req(v$active_pk)
      chr_offset <- d$chr_map[Chr == v$active_pk$Chr, Offset_Mbp]
      local_min <- max(0, x_range[1] - chr_offset)
      local_max <- min(d$chr_map[Chr == v$active_pk$Chr, Length / 1e6], x_range[2] - chr_offset)
      local_breaks <- pretty(c(local_min, local_max), n = 6)
      local_breaks <- local_breaks[local_breaks >= local_min & local_breaks <= local_max]
      global_breaks <- local_breaks + chr_offset
      diff_mbp <- local_max - local_min
      dec <- if (diff_mbp <= 0.05) 3 else if (diff_mbp <= 0.5) 2 else 1
      local_labels <- sprintf(paste0("%.", dec, "f"), local_breaks)
      p <- p + 
        scale_x_continuous(limits=x_range, breaks=global_breaks, labels=local_labels, expand=c(0,0)) +
        do.call(scale_y_continuous, y_scale_args) +
        labs(x = paste0(v$active_pk$Chr, " (Mbp)"))
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
             showlegend = !isTRUE(input$hide_manhattan_legend),
             margin = list(t = 50),
             xaxis = list(ticks = "outside", ticklen = 5, tickcolor = "black"),
             yaxis = list(ticks = "outside", ticklen = 5, tickcolor = "black")) %>% 
      config(displayModeBar = FALSE)
  })

  output$schematic <- renderPlot({
    if (is.null(v$active_pk)) {
      return(
        ggplot() + 
          annotate("text", x = 0.5, y = 0.5, 
                   label = "Click any peak or coding SNP on the Manhattan plot above\nor search a gene in the top bar to inspect locus architecture", 
                   size = 5.2, color = "#6c757d", fontface = "italic") + 
          theme_void() + 
          theme(panel.background = element_rect(fill = "#fdfdfe", color = "#e9ecef", linewidth = 1))
      )
    }
    pk <- v$active_pk
    win <- get_win_kb() * 1000
    
    window_genes <- d$genes[Chr == pk$Chr & End >= (pk$Mid - win) & Start <= (pk$Mid + win)]
    if(nrow(window_genes) == 0) return(NULL)
    
    plot_df <- merge(window_genes, d$degs[, .(GeneId, Category)], by="GeneId", all.x=TRUE)
    plot_df[is.na(Category), Category := "Non-DE"]
    if (!is.null(input$show_cat)) {
      plot_df[!Category %in% c(input$show_cat, "Non-DE"), Category := "Non-DE"]
    }
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
      annotate("text", x = pk$Mid, y = 1.55, label = paste("MafA peak", sub(".*peak_", "", pk$PeakID)), 
               color = "gold4", fontface = "bold", size = 6) +
      geom_text(aes(x = (Start + End)/2, y = y_lev + rect_h + 0.08, label = Symbol, color = Category), fontface = "bold", size = 6) +
      scale_fill_manual(name = "Category / Impact", values=all_fill_colors) + 
      scale_color_manual(name = "Category / Impact", values=all_color_colors) +
      scale_x_continuous(limits = c(pk$Mid - win, pk$Mid + win), labels = function(x) format(x/1e6, digits=5), expand = c(0, 0)) +
      labs(x=paste(pk$Chr, "(Mbp)"), y="") +
      theme_minimal() + 
      theme(
        axis.text.y=element_blank(), 
        panel.grid=element_blank(), 
        axis.line.x = element_line(),
        axis.text.x = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        legend.position = if (isTRUE(input$hide_schematic_legend)) "none" else "right"
      )
  })
  
  output$metadata_panel <- renderUI({
    req(v$active_pk); pk <- v$active_pk
    win_size <- get_win_kb() * 1000
    local_degs <- d$degs[Chr == pk$Chr & abs(Start - pk$Mid) <= win_size]
    if (!is.null(input$show_cat)) {
      local_degs <- local_degs[Category %in% input$show_cat]
    }
    
    selected_impacts <- if (is.null(input$snp_impact)) c("HIGH", "MODERATE") else input$snp_impact
    min_pc <- if (is.null(input$min_phastcons)) 0 else input$min_phastcons
    
    local_snps <- d$snps[Chr == pk$Chr & abs(pos - pk$Mid) <= win_size &
                           impact_simple %in% selected_impacts & phastCons_score >= min_pc]
    
    wellPanel(class = "well-meta",
      div(style = "display: flex; justify-content: space-between; align-items: center; border-bottom: 1.5px solid #e0c868; padding-bottom: 10px; margin-bottom: 14px;",
        div(class = "meta-title", style = "margin-bottom: 0;",
            paste("Detailed Peak Information:", sub(".*peak_", "", pk$PeakID))),
        div(style = "display: flex; gap: 18px; align-items: center;",
          span(strong("Location: "), sprintf("%.3f Mbp (%s)", pk$Mid/1e6, pk$Chr)),
          span(strong("Strain Divergent SNPs: "), pk$variant_count),
          if (!is.null(v$active_qtl) && v$active_qtl$Chr == pk$Chr && 
              pk$Mid >= v$active_qtl$ci.low * 1e6 && pk$Mid <= v$active_qtl$ci.high * 1e6) {
            span(style = "background-color: #6f42c1; color: white; padding: 4px 10px; border-radius: 4px; font-weight: bold; font-size: 0.88em;",
                 paste("Within F2 QTL:", v$active_qtl$trait, paste0("(LOD ", round(v$active_qtl$lod, 1), ")")))
          }
        )
      ),
      div(style = "display: flex; gap: 25px;",
        div(style = "flex: 1; border-right: 1px solid #eee; padding-right: 20px;",
          p(strong("Locus Differentially Expressed Genes (DEGs & log2FC):"), style = "margin-bottom: 8px; color: #333;"),
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
          } else { p("No DEGs in window.", style="font-style:italic; color: #777;") }
        ),
        div(style = "flex: 1;",
          p(strong("Coding SNPs in Locus:"), style = "margin-bottom: 8px; color: #333;"),
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
          } else { p("No coding SNPs in window.", style="font-style:italic; color: #777;") }
        )
      )
    )
  })
}

shinyApp(ui, server)
