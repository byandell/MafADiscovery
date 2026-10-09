# Developer Guide — MafA Discovery

Developer guide, architectural reference, and coding standards for the **MafA Discovery Integrated Genomic Explorer** repository.

* [Repository Architecture](#repository-architecture)
* [Data Layer & Coordinate System](#data-layer--coordinate-system)
* [Application Architecture (`SHINY_APP/app.R`)](#application-architecture-shiny_appappr)
  * [Data Engine (`prepare_data()`)](#data-engine-prepare_data)
  * [Reactive State Management](#reactive-state-management)
  * [Plotly Manhattan Plot & Dual-Axis Scaling](#plotly-manhattan-plot--dual-axis-scaling)
  * [Locus Schematic Rendering](#locus-schematic-rendering)
  * [Metadata Summary Panel](#metadata-summary-panel)
* [Shinylive & Serverless Deployment](#shinylive--serverless-deployment)
* [Coding Guidelines & Safety Rules](#coding-guidelines--safety-rules)
* [Version Control Governance](#version-control-governance)

---

## Repository Architecture

```
MafADiscovery/
├── app.R                          # Top-level launcher: shiny::runApp("SHINY_APP")
├── SHINY_APP/                     # Self-contained Shiny application source
│   ├── app.R                      # Main application UI and Server logic
│   ├── interpretation_guide.md    # In-app interpretation modal content
│   ├── developer_guide.md         # In-app developer guide modal content with links
│   ├── MafA_Peaks_with_SNPs_v3.csv              # MafA ChIP/CUT&RUN peaks and scores
│   ├── Master_DEG_Strain_Comparison_v3.csv      # Differential gene expression dataset
│   ├── mouse_genes_mm39_v3.csv                  # Ensembl GRCm39 gene backbone
│   ├── B6_SJL_prioritized_protein_coding_SNPs.csv # Coding SNPs with phastCons scores
│   ├── Top_glycemic_QTL_for_sex_additive_analysis.csv # F2 glycemic QTL loci with 95% CIs
│   └── README.md                  # Detailed collaborator documentation
├── .github/
│   └── workflows/
│       └── deploy-shinylive.yaml  # GitHub Actions automated Shinylive export & deployment
├── docs/                          # Standalone HTML documentation pages for GitHub Pages
│   ├── DEVELOPER.html             # Rendered master Developer Guide
│   ├── shinyapp.html              # Rendered Legacy Prototypes guide
│   ├── publishapp.html            # Rendered Publishing & WebAssembly guide
│   ├── redesign.html              # Rendered UI Redesign guide
│   ├── qtlanalysis.html           # Rendered QTL Analysis guide
│   └── .nojekyll                  # Bypasses Jekyll on GitHub Pages
├── render_docs.R                  # Generator script rendering *.md to docs/*.html
├── MafADiscovery.Rproj            # RStudio project configuration
├── shinyapp.md                    # Root architectural module: legacy prototypes
├── publishapp.md                  # Root architectural module: publishing guide
├── redesign.md                    # Root architectural module: UI redesign
├── qtlanalysis.md                 # Root architectural module: QTL integration
├── DEVELOPER.md                   # Root master developer guide (this file)
├── AGENTS.md                      # AI assistant project guidelines
├── LICENSE                        # MIT License
└── README.md                      # Project overview & quick start
```

---

## Data Layer & Coordinate System

All genomic coordinates in this project are standardized on the mouse **GRCm39 / mm39** genome assembly:

1. **Chromosome Mapping (`chr_map`)**:
   - Chromosomes: `Chr1` through `Chr19`, `ChrX`, and `ChrY`.
   - Lengths are standard GRCm39 base-pair lengths.
   - Global offsets (`Offset_Mbp` and `Midpoint_Mbp`) are computed via cumulative summation to place all 21 chromosomes along a single linear coordinate axis for the genome-wide Manhattan plot.
2. **MafA Binding Peaks (`MafA_Peaks_with_SNPs_v3.csv`)**:
   - Contains peak intervals (`Start`, `End`, `Mid`), peak scores, and variant counts (`variant_count`).
3. **Differential Gene Expression (`Master_DEG_Strain_Comparison_v3.csv`)**:
   - Classifies genes into 4 categories: `Shared`, `C57_Specific`, `SJL_Specific`, and `Discordant`.
   - Includes directional regulation (`UP` or `DOWN`) and strain-specific log2 fold-changes (`log2FC_C57`, `log2FC_SJL`).
4. **Gene Backbone (`mouse_genes_mm39_v3.csv`)**:
   - Comprehensive Ensembl GRCm39 annotation containing `GeneId`, `Symbol`, coordinates, and `Strand` (`1` for positive, `-1` for negative).
5. **Prioritized Coding SNPs (`B6_SJL_prioritized_protein_coding_SNPs.csv`)**:
   - High- and moderate-impact coding variants between C57BL/6J and SJL/J.
   - Contains evolutionary conservation scores (`phastCons_score`), variant consequences (`csq`), and amino acid alterations (`aa_change`).
6. **F2 Glycemic QTL Loci (`Top_glycemic_QTL_for_sex_additive_analysis.csv`)**:
   - Top glycemic QTLs from the B6 x SJL F2 study across traits, markers, chromosomes, peak positions (`pos`), 95% confidence intervals (`ci.low` to `ci.high`), and additive effects (`BB`, `BS`, `SS`).

---

## Application Architecture (`SHINY_APP/app.R`)

### Data Engine (`prepare_data()`)

The data engine runs once at session initialization:
- Loads the 5 CSV files using fast `data.table::fread()`.
- Calculates linear Manhattan positions (`GlobalPos_Mbp`) for peaks, coding variants, and QTL confidence intervals using `chr_map`.
- Formats SNP impacts (`HIGH` vs `MODERATE`) and ensures missing `phastCons_score` values default to `0`.
- Eliminates heavy Bioconductor packages (`GenomicRanges`, `IRanges`, `S4Vectors`) by using `data.table` rolling joins (`roll = "nearest"`):
  ```r
  res <- deg_cols[pks_q, on = .(Chr, deg_start = query_pos), roll = "nearest"]
  res <- res[abs(Mid - actual_deg_start) <= (input$win_kb * 1000) & !is.na(Category)]
  ```
  This reduces WebAssembly bundle sizes and dramatically speeds up browser startup.

### Reactive State Management

State is managed via `reactiveValues` in `v`:
- `v$active_pk`: Currently selected MafA binding peak row.
- `v$active_snp`: Currently selected coding SNP row (if clicked).
- `v$active_qtl`: Currently selected F2 glycemic QTL row.
- `v$last_gene`: Cached symbol of the last searched gene.
- `v$current_chr`: Active chromosome tracking (synced when peaks, genes, or QTLs are selected to allow seamless transitions between view modes).
- `v$user_zoom`: User-defined x-axis range captured from Plotly `relayout` events.
- `v$reset_trigger`: Integer incremented to reset Plotly view revision state.

### UI & Navigation Lifecycle

- **Browser Title & Metadata**: `fluidPage(title = "MafA Discovery: Integrated Genomic Explorer")` with `<title>` tag and inline DNA SVG favicon (`🧬`) in `tags$head`.
- **Top Global Navigation Bar**: A unified header card containing:
  - Brand header (`🧬 MafA Discovery: Genomic Explorer`).
  - Gene symbol search (`selectizeInput("search_gene", ...)`).
  - View mode dropdown (`selectInput("zoom_mode", ...)`).
  - Contextual selector container (`uiOutput("context_selector_ui")`): dynamically renders chromosome selector (`sel_chr`) in *Chromosome* mode, QTL selector (`sel_qtl`) in *QTL Region* mode, or locus indicators.
  - Quick action buttons: `📊 QTLs` modal, `ℹ️ Guide` interpretation modal, direct `📖 Dev Guide ↗` external link, and `↺ Reset`.
- **Contextual Plot Controls**:
  - **Manhattan Plot**: Adjacent control panel hosting legend visibility toggle (`hide_manhattan_legend`), visible peak category checkboxes (`show_cat`), min peak score filter (`min_score`), coding SNPs toggle (`show_snps_main`), and coding variant impact/phastCons filters (`snp_impact`, `min_phastcons`).
  - **Locus Schematic**: Adjacent control panel hosting legend visibility toggle (`hide_schematic_legend`), locus window size dropdown (`win_kb`), and a micro-architecture visual interpretation legend.
- **Genome-Wide Auto-Reset**: Selecting `"Genome-Wide"` in `zoom_mode` or clicking `↺ Reset` invokes a centralized `reset_to_default_state()` routine, clearing active peaks, SNPs, searched genes, user zooms, and returning all filters and options to their pristine default settings.
- **Legend Visibility Controls**: Discrete checkboxes (`hide_manhattan_legend`, `hide_schematic_legend`) dynamically suppress legends via Plotly `layout(showlegend = ...)` and ggplot `theme(legend.position = ...)`.
- **QTL Table Reference**: A dedicated modal (`input$show_qtl_table`) allows inspecting all 11 F2 glycemic QTLs with their LOD scores and additive effect estimates, and jumping directly to any QTL's confidence interval.
- **Interpretation Guide Modal**: Loaded dynamically from [`SHINY_APP/interpretation_guide.md`](SHINY_APP/interpretation_guide.md) via `render_markdown_file()`. Explains biological background, macro/micro panel coordination, coding variant impact levels, and navigation modes, allowing text updates without modifying R code.
- **Documentation Hub**: The `📖 Dev Guide ↗` button connects directly to the rendered documentation portal at `https://byandell.github.io/MafADiscovery/docs` (backed by `docs/index.html` compiled from [`SHINY_APP/developer_guide.md`](SHINY_APP/developer_guide.md)), linking to all modules and [`SHINY_APP/about.md`](SHINY_APP/about.md).

### Plotly Manhattan Plot & Dual-Axis Scaling

- Renders genome-wide, chromosome-level, QTL-level, or locus-level peak points using `ggplot2` and converts via `plotly::ggplotly()`.
- Automatically highlights the active F2 QTL confidence interval (`ci.low` to `ci.high`) with a shaded rectangular region and vertical dashed line at the QTL peak position.
- Dynamically scales local Mbp tick marks (`scale_x_continuous`) in *QTL Region* and *Locus Zoom* modes based on the active interval.
- Supports secondary y-axis scaling for coding SNPs:
  - Peak score plotted on primary y-axis.
  - Coding SNPs scaled dynamically to `phastCons_score * y_max` with `sec_axis(~ . / y_max, name = "phastCons Conservation Score (0-1)")`.
- Clicking either a peak or a coding variant triggers `plotly_click` handling, selecting the entity, highlighting it with a gold diamond, and opening *Locus Zoom*.
- Supports discrete legend visibility toggle (`hide_manhattan_legend`).

### Locus Schematic Rendering

When an active peak is selected (`v$active_pk`):
- Gene models within `± win_kb` are extracted and assigned vertical row offsets (`y_rows`) to avoid visual collisions.
- Transcription start sites (TSS) are computed based on strand orientation (`Strand == 1` vs `Strand == -1`).
- Directional arrow segments (`geom_segment(..., arrow = ...)`) indicate gene transcription direction.
- Coding variants are mapped to their respective gene models and rendered with distinct symbols and impact colors (`#CC00CC` for HIGH, `#DAA520` for MODERATE).
- Shows an informative guidance placeholder when no peak is selected.
- Supports discrete legend visibility toggle (`hide_schematic_legend`).

### Detailed Peak Information Card

Rendered below the Locus Schematic on the main canvas across the full 12-column grid via `shiny::renderUI()`, utilizing a wide two-column layout:
- **Header**: Active peak ID, Mbp coordinate on chromosome, variant count, and F2 QTL association badge.
- **Left Column**: Local DEGs with individual strain log2 fold-changes and direction.
- **Right Column**: Local coding variants grouped by gene, summarizing HIGH/MODERATE counts, consequence annotations, and `phastCons` scores.

---

## Modular Architectural Documentation

The repository maintains four specialized architectural modules at root, automatically synchronized to GitHub Pages during CI deployment:

1. **[`shinyapp.md`](shinyapp.md) — Legacy Standalone Prototypes**:
   - Comprehensive technical specifications and lineage for Version 1 (`MafA_Discovery_App.R`) and Version 2 (`MafA_Discovery_App_v2.R`).
   - Documents the original Bioconductor `GenomicRanges` implementation and data requirements.
2. **[`publishapp.md`](publishapp.md) — Publishing & Deployment**:
   - Architecture for static WebAssembly distribution via Shinylive (webR).
   - Automated GitHub Actions deployment (`deploy-shinylive.yaml`) and GitHub Pages hosting configuration.
3. **[`redesign.md`](redesign.md) — UI Redesign & Reactive Lifecycle**:
   - Dynamic view mode transitions (`Genome-Wide`, `Chromosome`, `QTL Region`, `Locus Zoom`).
   - Conditional control rendering, coordinate auto-scaling, and state persistence rules.
4. **[`qtlanalysis.md`](qtlanalysis.md) — F2 Glycemic QTL Integration**:
   - Ingestion of F2 study glycemic loci (`Top_glycemic_QTL_for_sex_additive_analysis.csv`).
   - Confidence interval auto-bounding, visual interval banners, and interactive QTL reference table.

---

## Shinylive & Serverless Deployment

The app is deployed to GitHub Pages as a static WebAssembly bundle:

1. **Workflow (`.github/workflows/deploy-shinylive.yaml`)**:
   - Triggers on push to `main` / `master` when files in `SHINY_APP/**`, root markdown guides (`*.md`), or the workflow itself change.
   - Sets up R on Ubuntu, installs `shiny`, `shinylive`, `data.table`, `ggplot2`, and `plotly`.
   - Executes `shinylive::export(appdir = "SHINY_APP", destdir = "site")`.
   - Copies root documentation assets (`DEVELOPER.md`, `shinyapp.md`, `publishapp.md`, `redesign.md`, `qtlanalysis.md`) into `site/` for public hosting.
   - Uploads `site/` and deploys to GitHub Pages via `actions/deploy-pages@v4`.
2. **Local Static Verification**:
   - If testing a static export locally with WebAssembly:
     ```r
     shinylive::export(appdir = "SHINY_APP", destdir = "site")
     httpuv::runStaticServer("site", port = 8888)
     ```

---

## Coding Guidelines & Safety Rules

- **Vector Subsetting Safety**: Always use `grepl("^\\s*#'", lines)` with `!grepl(...)` or `grep(..., invert = TRUE)`. Never use `!grep(...)` in R.
- **Explicit Package Namespacing**: Use explicit namespaces (`data.table::fread`, `shiny::renderPlot`, `plotly::ggplotly`) where appropriate.
- **Portability**: All data file paths inside `app.R` must remain relative so that the app executes identically in local R sessions, Shinylive WebAssembly, and Docker/cloud deployments.

---

## Version Control Governance

- **No Automatic Git Commit/Push**: Prepare file edits and run local checks, but leave staging, committing, and pushing for manual user execution.
- **Preserve Data Files**: Never overwrite or delete raw CSV datasets destructively.
