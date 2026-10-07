# MafA Discovery: Integrated Genomic Explorer

An interactive R Shiny application for exploring MafA transcription factor binding peaks, strain-divergent genetic variants between C57BL/6J (B6) and SJL/J mice, and downstream differential gene expression (DEGs).

* [Overview](#overview)
* [Repository Structure](#repository-structure)
* [Prerequisites & Installation](#prerequisites--installation)
* [Launching the Application](#launching-the-application)
* [Static Web Deployment (GitHub Pages & Shinylive)](#static-web-deployment-github-pages--shinylive)
* [Data Layer Catalog](#data-layer-catalog)
* [User Interface & Key Features](#user-interface--key-features)
* [Developer Documentation](#developer-documentation)
* [License](#license)

---

## Overview

The **MafA Discovery Explorer** integrates multiple layers of genomic and molecular data from C57BL/6J (B6) and SJL/J mouse models:

1. **MafA Binding Peaks**: ChIP-seq / CUT&RUN binding peak locations in GRCm39 coordinates with binding intensity scores.
2. **Strain-Divergent Variants (SNPs & Indels)**: Variant counts between C57BL/6J and SJL/J localized within MafA binding intervals.
3. **Differential Gene Expression (DEGs)**: Strain comparisons between C57BL/6J and SJL/J with log2 fold-changes, directional shifts, and categorical classifications (`Shared`, `C57_Specific`, `SJL_Specific`, `Discordant`).
4. **Prioritized Protein-Coding Variants**: High- and moderate-impact coding sequence polymorphisms between B6 and SJL with evolutionary conservation scores (`phastCons`), consequences, and amino acid substitutions.

This tool allows researchers and collaborators to rapidly scan genome-wide associations, examine local chromatin architecture around candidate genes, and evaluate whether sequence polymorphisms in MafA binding regions or coding sequences correlate with strain-divergent gene expression.

---

## Repository Structure

```
MafADiscovery/
├── app.R                          # Root launcher (runs shiny::runApp("SHINY_APP"))
├── SHINY_APP/                     # Production Shiny application directory
│   ├── app.R                      # Main self-contained Shiny application script
│   ├── MafA_Peaks_with_SNPs_v3.csv              # MafA binding peak calls & scores
│   ├── Master_DEG_Strain_Comparison_v3.csv      # Differential gene expression dataset
│   ├── mouse_genes_mm39_v3.csv                  # Ensembl GRCm39 gene annotation backbone
│   ├── B6_SJL_prioritized_protein_coding_SNPs.csv # Coding SNPs with phastCons scores
│   └── README.md                  # Detailed collaborator user guide
├── .github/
│   └── workflows/
│       └── deploy-shinylive.yaml  # Automated CI/CD Shinylive export & Pages deployment
├── docs/                          # Shinylive static web distribution bundle
├── MafADiscovery.Rproj            # RStudio project configuration
├── shinyapp.md                    # Shinylive & GitHub Pages technical notes
├── DEVELOPER.md                   # Architecture & developer reference guide
├── AGENTS.md                      # AI assistant project guidelines
├── LICENSE                        # MIT License
└── README.md                      # Repository overview (this file)
```

---

## Prerequisites & Installation

The application requires **R** (>= 4.0) and relies exclusively on standard CRAN packages:

```r
install.packages(c("shiny", "data.table", "ggplot2", "plotly"))
```

*Note: The application uses optimized `data.table` rolling joins instead of heavy Bioconductor packages, allowing it to load quickly and run serverlessly in web browsers via WebAssembly/Shinylive.*

---

## Launching the Application

### Option A: From RStudio (Recommended)

1. Open [`MafADiscovery.Rproj`](MafADiscovery.Rproj) in RStudio.
2. Open [`app.R`](app.R) or [`SHINY_APP/app.R`](SHINY_APP/app.R).
3. Click the **Run App** button in the top-right corner of the editor toolbar, or run:
   ```r
   shiny::runApp()
   ```

### Option B: From the R Console or Terminal

Launch directly from the repository root:

```r
shiny::runApp("SHINY_APP")
```

Or run from the command line:

```bash
R -e 'shiny::runApp("SHINY_APP")'
```

---

## Static Web Deployment (GitHub Pages & Shinylive)

This application can run 100% serverlessly in any modern web browser using **Shinylive** (webR):

1. **Automated CI/CD Workflow**: [`.github/workflows/deploy-shinylive.yaml`](.github/workflows/deploy-shinylive.yaml) automatically builds and exports the app using `shinylive::export()` on every push to `main` or `master`.
2. **GitHub Pages Configuration**:
   - In GitHub repository settings: **Settings** $\to$ **Pages**.
   - Under **Build and deployment** > **Source**, choose **GitHub Actions**.
   - The live application will be published automatically at `https://<username>.github.io/MafADiscovery/`.

---

## Data Layer Catalog

The application relies on four curated CSV datasets stored in [`SHINY_APP/`](SHINY_APP/):

| File | Description | Key Fields |
| :--- | :--- | :--- |
| [`MafA_Peaks_with_SNPs_v3.csv`](SHINY_APP/MafA_Peaks_with_SNPs_v3.csv) | MafA ChIP/CUT&RUN peak intervals and intensity scores | `PeakID`, `Chr`, `Start`, `End`, `Mid`, `Peak Score`, `variant_count`, `Dist_TSS_Num` |
| [`Master_DEG_Strain_Comparison_v3.csv`](SHINY_APP/Master_DEG_Strain_Comparison_v3.csv) | B6 vs SJL differential expression statistics | `GeneId`, `Symbol`, `Chr`, `Start`, `End`, `Category`, `Direction`, `log2FC_C57`, `log2FC_SJL` |
| [`mouse_genes_mm39_v3.csv`](SHINY_APP/mouse_genes_mm39_v3.csv) | Ensembl GRCm39 mouse gene genomic coordinates | `GeneId`, `Symbol`, `Chr`, `Start`, `End`, `Strand`, `Type` |
| [`B6_SJL_prioritized_protein_coding_SNPs.csv`](SHINY_APP/B6_SJL_prioritized_protein_coding_SNPs.csv) | High- and moderate-impact coding sequence SNPs | `variant_id`, `chr`, `pos`, `gene_symbol_1`, `csq`, `aa_change`, `phastCons_score`, `impact` |

---

## User Interface & Key Features

* **Genome-Wide Manhattan View**: Visualizes MafA peak scores across all 21 mouse chromosomes (1–19, X, Y), overlaid with coding sequence SNPs sized and shaded by `phastCons` conservation scores.
* **Interactive Navigation & Zoom**: Toggle seamlessly between *Genome-Wide*, *Chromosome*, and *Locus Zoom* modes, or search directly for any mouse gene symbol via the searchable autocomplete dropdown.
* **Locus Schematic View**: Detailed local genomic view centered on the active MafA binding peak. Displays gene structures with strand-directed transcription start site (TSS) arrows, nearby coding SNPs colored by predicted impact, and DEG status color coding.
* **Locus Metadata Panel**: Interactive summary card reporting peak coordinates, variant counts, log2 fold changes for all local DEGs, and detailed consequence and amino acid changes for local coding variants.

---

## Developer Documentation

For technical architectural details, data engine documentation, and coding guidelines, see [**`DEVELOPER.md`**](DEVELOPER.md).

---

## License

This project is licensed under the [MIT License](LICENSE).
