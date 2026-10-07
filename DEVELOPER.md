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
│   ├── MafA_Peaks_with_SNPs_v3.csv              # MafA ChIP/CUT&RUN peaks and scores
│   ├── Master_DEG_Strain_Comparison_v3.csv      # Differential gene expression dataset
│   ├── mouse_genes_mm39_v3.csv                  # Ensembl GRCm39 gene backbone
│   ├── B6_SJL_prioritized_protein_coding_SNPs.csv # Coding SNPs with phastCons scores
│   └── README.md                  # Detailed collaborator documentation
├── .github/
│   └── workflows/
│       └── deploy-shinylive.yaml  # GitHub Actions automated Shinylive export & deployment
├── docs/                          # Shinylive static web distribution bundle (gitignored)
├── MafADiscovery.Rproj            # RStudio project configuration
├── shinyapp.md                    # Shinylive deployment notes & configuration
├── DEVELOPER.md                   # Technical reference & developer guide (this file)
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

---

## Application Architecture (`SHINY_APP/app.R`)

### Data Engine (`prepare_data()`)

The data engine runs once at session initialization:
- Loads the 4 CSV files using fast `data.table::fread()`.
- Calculates linear Manhattan positions (`GlobalPos_Mbp`) for both peaks and coding variants using `chr_map`.
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
- `v$user_zoom`: User-defined x-axis range captured from Plotly `relayout` events.
- `v$reset_trigger`: Integer incremented to reset Plotly view revision state.

### Plotly Manhattan Plot & Dual-Axis Scaling

- Renders genome-wide or chromosome-level peak points using `ggplot2` and converts via `plotly::ggplotly()`.
- Supports secondary y-axis scaling for coding SNPs:
  - Peak score plotted on primary y-axis.
  - Coding SNPs scaled dynamically to `phastCons_score * y_max` with `sec_axis(~ . / y_max, name = "phastCons Conservation Score (0-1)")`.
- Clicking either a peak or a coding variant triggers `plotly_click` handling, selecting the entity and highlighting it with a gold diamond.

### Locus Schematic Rendering

When an active peak is selected (`v$active_pk`):
- Gene models within `± win_kb` are extracted and assigned vertical row offsets (`y_rows`) to avoid visual collisions.
- Transcription start sites (TSS) are computed based on strand orientation (`Strand == 1` vs `Strand == -1`).
- Directional arrow segments (`geom_segment(..., arrow = ...)`) indicate gene transcription direction.
- Coding variants are mapped to their respective gene models and rendered with distinct symbols and impact colors (`#CC00CC` for HIGH, `#DAA520` for MODERATE).

### Metadata Summary Panel

Renders HTML cards via `shiny::renderUI()` summarizing:
- Active peak coordinates, length, and variant count.
- Local DEGs with individual strain log2 fold-changes.
- Local coding variants with gene symbols, consequence annotations, and `phastCons` conservation scores.

---

## Shinylive & Serverless Deployment

The app is deployed to GitHub Pages as a static WebAssembly bundle:

1. **Workflow (`.github/workflows/deploy-shinylive.yaml`)**:
   - Triggers on push to `main` / `master` when files in `SHINY_APP/**` or the workflow itself change.
   - Sets up R on Ubuntu, installs `shiny`, `shinylive`, `data.table`, `ggplot2`, and `plotly`.
   - Executes `shinylive::export(appdir = "SHINY_APP", destdir = "site")`.
   - Uploads `site/` and deploys to GitHub Pages via `actions/deploy-pages@v4`.
2. **Local Static Verification**:
   - To test the static build locally with WebAssembly:
     ```r
     httpuv::runStaticServer("docs", port = 8888)
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
