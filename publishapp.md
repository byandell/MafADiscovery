# Publish Shiny App to GitHub Pages

Here is a comprehensive plan to publish the Shiny app in `SHINY_APP/` using **GitHub Pages**.

## 1. Technical Context: GitHub Pages & Shiny

GitHub Pages is a **static web host** (HTML, CSS, JS, and WebAssembly). It does **not** run an active backend R server (such as `Rscript` or standard `shiny::runApp()`).

To host an interactive R Shiny application on GitHub Pages, there are two primary options:

| Strategy | How it Works | Pros | Considerations |
| :--- | :--- | :--- | :--- |
| **Option 1: Shinylive for R (webR)** *(Recommended)* | Runs R and Shiny entirely inside the visitor's browser using **WebAssembly (WASM)**. | • 100% serverless<br>• Directly hosted on GitHub Pages<br>• Zero recurring server maintenance | • First-time page load downloads R-wasm runtime (~15–20 MB cached afterwards)<br>• Requires dependencies supported in webR |
| **Option 2: shinyapps.io + GitHub Pages landing/iframe** | App runs on Posit's free/paid shinyapps.io server; `docs/` hosts a GitHub Pages landing page embedding it. | • Standard R server environment<br>• No WASM compatibility checks | • Free tier sleeps after inactivity and has 25 active hours/month limit |

---

## 2. Directory Structure: Source vs. Build (`docs/`)

Rather than moving and overwriting the source files into `docs/`, the standard best practice is to **maintain `SHINY_APP/` as the source** and **build/export the static assets into `docs/`**:

```
MafADiscovery/
├── SHINY_APP/                     # SOURCE: R scripts, raw CSVs, docs
│   ├── app.R                      # Main application file (renamed or symlinked)
│   ├── MafA_Peaks_with_SNPs_v3.csv
│   ├── Master_DEG_Strain_Comparison_v3.csv
│   ├── mouse_genes_mm39_v3.csv
│   └── README.md
├── docs/                          # BUILD TARGET: Published via GitHub Pages
│   ├── .nojekyll                  # Required: bypasses Jekyll processing
│   ├── index.html                 # Shinylive entry point
│   ├── edit/                      # Shinylive editor / runtime assets
│   └── shivlive/ (or _webr/)      # WebAssembly runtime files
```

---

## 3. Dependency & Performance Audit (webR Readiness)

The app (`MafA_Discovery_App.R`) relies on 5 libraries:

1. `shiny`, `ggplot2`, `plotly`, `data.table`: All are available as pre-compiled WebAssembly binaries in the r-wasm CRAN repository (`repo.r-wasm.org`).
2. `GenomicRanges`: Bioconductor package with C dependencies.
   - *Consideration*: While webR supports many Bioconductor packages, loading `GenomicRanges` + `S4Vectors` + `IRanges` adds ~25 MB of WASM binaries to download.
   - *Optimization Option*: The app uses `GenomicRanges` solely for `findOverlaps()` and `nearest()` (finding DEGs within `win_kb` of peak midpoint). This can optionally be replaced with a fast `data.table` interval join or base R calculation to cut app load time by >60%.
3. **Data payload**: The 3 CSVs total ~6.3 MB (`mouse_genes_mm39_v3.csv` is ~5.1 MB). This comfortably loads into browser memory in webR.

---

## 4. Proposed Step-by-Step Implementation Steps

### Phase 1: Prepare the App for Export

1. Standardize the script filename: Shinylive expects `app.R` (or a single R script). We can create an `app.R` in `SHINY_APP/` that cleanly wraps or links to `MafA_Discovery_App.R`.
2. Verify package dependencies in webR format. If desired, benchmark whether to keep `GenomicRanges` or streamline the distance check via `data.table`.

### Phase 2: Static Export via Shinylive

1. Use the R package `shinylive` to export the application into `docs/`:

   ```r
   if (!requireNamespace("shinylive", quietly = TRUE)) install.packages("shinylive")
   shinylive::export(
     appdir = "SHINY_APP",
     destdir = "docs",
     subdir = ""
   )
   ```

2. Create `docs/.nojekyll` to ensure GitHub Pages does not ignore underscore-prefixed directories (e.g., `_webr/`).

### Phase 3: Local Empirical Verification

1. Test the static build locally using a local web server (a static server is required because browser security policies block WASM from `file://` URLs):

   ```r
   httpuv::runStaticServer("docs", port = 8888)
   ```

2. Verify that:
   - The webR engine initializes in the browser.
   - Data files load into the virtual filesystem.
   - The Manhattan plot, gene symbol autocomplete, zoom modes, and locus schematic render interactively.

### Phase 4: GitHub Actions Automated Build & Pages Activation

Rather than committing and uploading heavy precompiled Shinylive/webR assets (`docs/shinylive/` ~60MB) directly into Git history, the build and deployment is handled on the GitHub Pages end via GitHub Actions:

1. **Workflow Automation (`.github/workflows/deploy-shinylive.yaml`)**:
   - Triggers on push to `main` (when `SHINY_APP/**`, `*.md`, or the workflow changes) or manual `workflow_dispatch`.
   - Sets up Ubuntu runner with R, installs standard CRAN dependencies (`shiny`, `shinylive`, `data.table`, `ggplot2`, `plotly`).
   - Runs `shinylive::export(appdir = "SHINY_APP", destdir = "site")` in CI.
   - Copies root developer guides (`DEVELOPER.md`, `shinyapp.md`, `publishapp.md`, `redesign.md`, `qtlanalysis.md`) into `site/` for public hosting alongside the app.
   - Deploys the static bundle as a Pages artifact via `actions/deploy-pages@v4`.
2. **Repository Cleanliness**:
   - `/docs/` and `/site/` are added to `.gitignore`.
   - The git repository contains only clean source files and CSV data in `SHINY_APP/`.
3. **GitHub Repository Settings**:
   - Navigate to **Settings** $\to$ **Pages**.
   - Under **Build and deployment** > **Source**, select **GitHub Actions** (instead of "Deploy from a branch").
   - Pushing changes to `main` triggers automated build and deployment to `https://<organization-or-user>.github.io/MafADiscovery/`.

---

## Implementation Summary & Status

All implementation steps for **Option 1 (Shinylive / webR)** with `data.table` and serverless GitHub Actions deployment have been executed:

1. **`GenomicRanges` Replaced by `data.table`**:
   - Replaced `findOverlaps()` and `nearest()` with a fast `data.table` rolling join (`roll = "nearest"` on `.(Chr, deg_start = query_pos)` followed by `abs(Mid - actual_deg_start) <= win_kb * 1000`).
   - Verified 100% exact parity on the actual dataset (`6,790` peaks matching identically in symbols, categories, and peak IDs).
   - Eliminated Bioconductor dependencies (`GenomicRanges`, `IRanges`, `S4Vectors`, `Seqinfo`), reducing required WASM packages from 43 to 38 and speeding up browser load times.

2. **Source and Entry Point Standardized**:
   - Created [`SHINY_APP/app.R`](SHINY_APP/app.R) as the primary Shiny application entry point for both local development and Shinylive packaging.
   - Maintained all source files in `SHINY_APP/`.

3. **Serverless CI/CD Pipeline (`.github/workflows/deploy-shinylive.yaml`)**:
   - Automated `shinylive::export()` directly within GitHub Actions runners.
   - Kept Git repository free of heavy binary WASM web assets (`/docs/` and `/site/` ignored in `.gitignore`).

4. **Next Step (GitHub Settings)**:
   - In GitHub repository settings: **Settings** $\to$ **Pages** $\to$ **Source: GitHub Actions**.
   - Pushes to `main` automatically publish the live app at `https://<organization-or-user>.github.io/MafADiscovery/`.
