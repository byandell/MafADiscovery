# MafA Discovery: Integrated Genomic Explorer

An interactive R Shiny application for exploring MafA transcription factor binding peaks, strain-divergent genetic variants between C57BL/6J and SJL/J mice, and downstream differential gene expression (DEGs).

* [Overview](#overview)
* [Application Architecture & `app.R`](#application-architecture--appr)
* [Prerequisites & Installation](#prerequisites--installation)
* [Directory & Data Files](#directory--data-files)
* [Launching the Application](#launching-the-application)
* [Static Web Deployment (GitHub Pages & Shinylive)](#static-web-deployment-github-pages--shinylive)
* [User Interface & Features](#user-interface--features)
  * [Genome-Wide Manhattan View](#genome-wide-manhattan-view)
  * [Navigation & Search](#navigation--search)
  * [Locus Schematic View](#locus-schematic-view)
  * [Filter Controls & Metadata Panel](#filter-controls--metadata-panel)
* [Biological Context & Interpretation](#biological-context--interpretation)
* [Legacy Standalone Applications](#legacy-standalone-applications)
* [Troubleshooting](#troubleshooting)

---

## Overview

The **MafA Discovery Explorer** integrates multiple layers of genomic and molecular data from the Vanderbilt C57BL/6J (B6) and SJL/J mouse models:

1. **MafA Binding Peaks**: ChIP-seq / CUT&RUN binding peak locations (GRCm39 coordinates) with binding peak intensity scores.
2. **Genetic Polymorphisms (SNPs & Indels)**: Variant counts between C57BL/6J and SJL/J localized within MafA binding intervals.
3. **Differential Gene Expression (DEGs)**: Strain comparisons between C57BL/6J and SJL/J with log2 fold-changes, directional shifts, and categorical classifications.
4. **Protein-Coding Variants**: High- and moderate-impact coding sequence polymorphisms between B6 and SJL with evolutionary conservation scores (`phastCons`), consequences, and amino acid substitutions.

This tool allows collaborators to rapidly scan genome-wide associations, examine local chromatin architecture around candidate genes, and evaluate whether sequence polymorphisms in MafA binding regions or coding sequences correlate with strain-divergent gene expression.

---

## Application Architecture & `app.R`

[`app.R`](app.R) serves as the primary, self-contained application entry point for both local execution and static Shinylive/GitHub Pages web deployment:

- **Integrated Feature Set**: Combines MafA binding peak analysis, DEG classification, and strain-divergent protein-coding SNPs with evolutionary conservation (`phastCons`) scoring.
- **Fast & Dependency-Free**: Replaces heavy Bioconductor dependencies (`GenomicRanges`, `IRanges`, `S4Vectors`) with optimized `data.table` rolling joins (`roll = "nearest"`). This yields 100% mathematical parity with earlier prototypes while enabling rapid in-browser execution via WebAssembly without heavy package downloads.
- **Standardized Entry Point**: Uses the standard Shiny file structure (`app.R`) for automated toolchains (`shinylive::export`), Shiny Server, Posit Connect, and one-click RStudio **Run App** execution.

> [!NOTE]
> Earlier standalone prototype scripts (`MafA_Discovery_App.R` and `MafA_Discovery_App_v2.R`) were historically developed in the Vanderbilt analysis repository. The current `app.R` consolidates all features with zero Bioconductor runtime dependencies for fast browser execution.

---

## Prerequisites & Installation

The application requires **R** (>= 4.0) and relies exclusively on standard CRAN packages:

```r
install.packages(c("shiny", "data.table", "ggplot2", "plotly"))
```

---

## Directory & Data Files

To run the application locally or export it, ensure the files reside together in the `SHINY_APP/` folder:

| File | Description |
| :--- | :--- |
| [`app.R`](app.R) | Primary production R Shiny application script (used for local launch and Shinylive export) |
| `MafA_Peaks_with_SNPs_v3.csv` | MafA binding peaks (`PeakID`, `Chr`, `Start`, `End`, `Mid`, `Peak Score`, `variant_count`, `Dist_TSS_Num`) |
| `Master_DEG_Strain_Comparison_v3.csv` | Differential expression data (`GeneId`, `Symbol`, `Chr`, `Start`, `End`, `Category`, `Direction`, `log2FC_C57`, `log2FC_SJL`) |
| `mouse_genes_mm39_v3.csv` | Mouse GRCm39 gene annotation backbone (`GeneId`, `Symbol`, `Chr`, `Start`, `End`, `Strand`, `Type`) |
| `B6_SJL_prioritized_protein_coding_SNPs.csv` | Prioritized coding SNPs with consequence, amino acid changes, and phastCons scores |
| `QTLresults/Top_glycemic_QTL_for_sex_additive_analysis.csv` | F2 glycemic QTL loci with 95% confidence intervals and additive effects |
| `QTLresults/scan_chr*.png` | Additive sex QTL LOD scans for AUC and trajectory glycemic traits |
| `interpretation_guide.md` | In-app modal documentation explaining macro/micro panels and navigation |
| `developer_guide.md` | In-app modal documentation linking to architectural guides and GitHub repo |
| [`README.md`](README.md) | Documentation for the active application |

> [!NOTE]
> The app loads the CSV files using relative file paths (e.g. `fread("Master_DEG_Strain_Comparison_v3.csv")`). Keep the CSV files in the same directory as `app.R`.

---

## Launching the Application

### Option A: From RStudio (Recommended)

1. Open `SHINY_APP/app.R` in RStudio.
2. Click the **Run App** button in the top-right corner of the editor toolbar, or run:
   ```r
   shiny::runApp()
   ```

### Option B: From the R Console

Set your working directory to the `SHINY_APP` folder and launch:

```r
setwd("SHINY_APP")
shiny::runApp()
```

Alternatively, launch directly from the repository root:

```r
shiny::runApp("SHINY_APP")
```

### Option C: From the Terminal

```bash
cd SHINY_APP
R -e "shiny::runApp()"
```

---

## Static Web Deployment (GitHub Pages & Shinylive)

The application is deployed to **GitHub Pages** via **Shinylive for R (webR)** using an automated **GitHub Actions** CI/CD workflow ([`.github/workflows/deploy-shinylive.yaml`](../.github/workflows/deploy-shinylive.yaml)):

1. **Automated CI/CD Build**:
   - Pushes to `main` affecting `SHINY_APP/` automatically trigger the GitHub Actions workflow.
   - The workflow spins up an Ubuntu runner, installs standard CRAN dependencies, and exports the static WebAssembly bundle via `shinylive::export(appdir = "SHINY_APP", destdir = "site")`.
   - The compiled site is uploaded as a GitHub Pages artifact and deployed directly via `actions/deploy-pages@v4`.
   - **Repository Cleanliness**: Precompiled WebAssembly assets (`shinylive/`, `app.json`) are never committed to Git; `/docs/` and `/site/` are listed in `.gitignore`.
2. **GitHub Pages Configuration**:
   - In repository **Settings** $\to$ **Pages**, set **Source: GitHub Actions**.
3. **Optional Local Static Testing**:
   If you wish to test static export locally without uploading artifacts:
   ```r
   shinylive::export(appdir = "SHINY_APP", destdir = "docs")
   httpuv::runStaticServer("docs")
   ```
   *(Local `docs/` is ignored by Git to avoid uploading vendor assets).*

---

## User Interface & Features

### Genome-Wide Manhattan View

The upper panel displays a Manhattan plot where each point represents a MafA binding peak mapped to nearby differentially expressed genes, with optional coding SNP overlays:

- **X-axis**: Genome position (Chr 1–19, Chr X, Chr Y) with cumulative offset spacing.
- **Left Y-axis**: MafA Peak Score (binding intensity/confidence).
- **Right Y-axis**: `phastCons` Conservation Score ($0–1$).
- **Point Fill Color (Peaks)**: Differential expression strain category:
  - 🟢 **Shared** (`#228B22`): Differentially expressed in both C57BL/6J and SJL/J in the same direction.
  - 🔵 **C57_Specific** (`#0000CC`): Differentially expressed exclusively in C57BL/6J.
  - 🔴 **SJL_Specific** (`#CC0000`): Differentially expressed exclusively in SJL/J.
  - 🟠 **Discordant** (`#FF8C00`): Differentially expressed in both strains but in opposite directions.
- **Point Shape (Peaks)**: Expression direction:
  - ▲ **UP**: Upregulated
  - ▼ **DOWN**: Downregulated
- **Point Size (Peaks)**: Number of polymorphic SNPs/indels within the peak interval between C57BL/6J and SJL/J.
- **Coding SNP Points**: Plotted as diamonds (`shape = 23`) scaled by conservation and colored by impact:
  - 🟣 **HIGH**: `#CC00CC` (stop-gain, frameshift, splice site)
  - 🟡 **MODERATE**: `#DAA520` (missense)
- **Active Selection**: The currently selected peak or coding SNP is highlighted with a large **gold diamond**.

### Navigation & Search

- **Search Gene Symbol**: Autocomplete search bar. Typing and selecting a gene symbol automatically finds the nearest MafA peak on that chromosome, selects it as active, and zooms into the locus.
- **View Mode Selector**:
  - **Genome-Wide**: Displays all 21 mouse chromosomes simultaneously.
  - **Chromosome**: Centers and scales on the chromosome selected in the dropdown.
  - **Locus Zoom**: Zooms into a window surrounding the active peak.
- **Interactive Plotly Controls**: Click directly on any peak point or coding SNP in the Manhattan plot to select it. Drag to zoom in on any region; double-click to reset.
- **↺ Reset to Genome-Wide**: Instantly clears active peak selection, resets manual zoom, and returns the view to the full genome.

### Locus Schematic View

When an active peak is selected, the lower panel renders a gene locus model centered on the peak:

- **Peak Interval**: Marked with a light yellow highlighted box (`#FFF9C4`) and a gold dashed vertical line at the peak midpoint.
- **Gene Bodies**: Plotted as horizontal bars across multiple non-overlapping vertical tiers.
  - Color-coded by DEG category (Shared, C57-specific, SJL-specific, Discordant).
  - Unperturbed protein-coding genes in the window are displayed in grey (`#D3D3D3`) for genomic context.
- **Transcription Start Sites (TSS) & Strand Direction**: For DEGs, directed arrows at the TSS indicate the strand of transcription (`+` forward or `-` reverse).
- **Coding SNP Markers**: Plotted directly on top of gene tracks:
  - Hollow diamonds (`shape = 5`) for lower conservation ($\text{phastCons} < 0.7$).
  - Solid diamonds (`shape = 23`) for high conservation ($\text{phastCons} \ge 0.7$).
- **Labels**: Gene symbols and the MafA peak ID are labeled directly above the features.

### Filter Controls & Metadata Panel

- **Visible Categories**: Checkboxes to toggle visibility of `Shared`, `C57_Specific`, `SJL_Specific`, and `Discordant` peaks.
- **Locus Window (Kbp)**: Numeric input (default `500` kb) defining the genomic radius used to map peaks to DEGs and determine the window size of the lower schematic plot.
- **Min Peak Score**: Slider to filter out peaks below a specified score threshold (range: 0 to 10,000).
- **Coding SNP Filters**:
  - Toggle SNP overlay on the main plot.
  - Impact checkboxes (`HIGH`, `MODERATE`).
  - Minimum `phastCons` slider ($0.0$ to $1.0$, default $0.7$).
- **Metadata Panel**: Displayed in the sidebar when an active peak or SNP is selected:
  - Peak identifier, chromosome position in Mbp, and variant count.
  - List of local DEGs with $\log_2\text{FC}$ values for C57 and SJL.
  - **Coding SNPs in Locus**: Breakdown of coding variants by gene symbol, highlighting high-conservation substitutions with exact positions, consequences, amino acid changes, and phastCons scores.

---

## Biological Context & Interpretation

- **MafA** is a master transcriptional regulator of pancreatic $\beta$-cell maturation and insulin gene expression.
- **Strain Divergence**: C57BL/6J and SJL/J strains show marked differences in $\beta$-cell function, glucose tolerance, and susceptibility to metabolic dysfunction.
- **Variant-Peak Overlap**: When a MafA peak contains a high count of strain-divergent SNPs (`variant_count > 0`), it may disrupt or alter transcription factor binding affinity, potentially driving strain-specific downstream gene regulation (`C57_Specific` or `SJL_Specific` DEGs).
- **Coding Variant Impact**: Coding sequence polymorphisms in candidate locus genes can cause functional amino acid alterations that synergize with regulatory variation to produce strain-specific physiological phenotypes.

---

## Historical Version Progression

Earlier iterations of the explorer were developed iteratively in the parent project:
- **v1 Prototype**: Baseline version focusing exclusively on MafA binding peaks and DEGs.
- **v2 Prototype**: Enhanced prototype introducing coding SNP overlays and dual Y-axis plotting using Bioconductor's `GenomicRanges`.
- **Production `app.R`**: Current unified architecture replacing Bioconductor with `data.table` rolling joins for fast WebAssembly execution, self-contained within this repository.

---

## Troubleshooting

- **`cannot open file '...csv': No such file or directory`**:
  Make sure your R working directory contains the CSV files, or launch the app using `shiny::runApp("SHINY_APP")` from the repository root.
- **Plotly zoom not updating**:
  Click the **↺ Reset to Genome-Wide** button in the sidebar to reset interactive Plotly coordinate states.
- **Local static web preview issues**:
  Because browser security settings disallow WebAssembly and service workers over raw `file://` URLs, always preview the static build using a web server such as `httpuv::runStaticServer("docs")`.
