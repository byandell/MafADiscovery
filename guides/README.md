# MafA Discovery — Documentation & Architecture Guides

Welcome to the `guides/` directory for the **MafA Discovery Integrated Genomic Explorer**. This folder consolidates all architectural documentation, in-app modal guides, technical implementation plans, and deployment specifications.

---

## 1. Guides Directory Catalog

| Guide | Target Audience | Description | Online HTML |
| :--- | :--- | :--- | :--- |
| [**`developer_guide.md`**](developer_guide.md) | Developers & Collaborators | Central index and documentation hub linking all architectural modules and GitHub source repositories. | [Documentation Hub](https://byandell.github.io/MafADiscovery/docs/index.html) |
| [**`user_guide.md`**](user_guide.md) | End Users & Researchers | In-app modal user guide explaining Manhattan macro view, Locus micro schematic, coding variant impacts, and navigation modes. | [User Guide](https://byandell.github.io/MafADiscovery/docs/user_guide.html) |
| [**`about.md`**](about.md) | All Users | Biological overview of MafA transcription factor dynamics, Vanderbilt/UW research consortium, and integrated data layers. | [About MafA Discovery](https://byandell.github.io/MafADiscovery/docs/about.html) |
| [**`qtlanalysis.md`**](qtlanalysis.md) | Geneticists & Developers | Technical specification of the F2 glycemic QTL integration, additive sex LOD scan plots, and interactive overview table. | [QTL Analysis](https://byandell.github.io/MafADiscovery/docs/qtlanalysis.html) |
| [**`redesign.md`**](redesign.md) | Developers & UI Engineers | Chronological implementation record of UI modernization, state reset lifecycle, contextual control placement, and mobile responsiveness. | [UI Redesign](https://byandell.github.io/MafADiscovery/docs/redesign.html) |
| [**`publishapp.md`**](publishapp.md) | DevOps & Developers | Static WebAssembly (WASM) Shinylive export strategy, dependency elimination, and GitHub Actions continuous deployment. | [Publishing & Deployment](https://byandell.github.io/MafADiscovery/docs/publishapp.html) |
| [**`shinyapp.md`**](shinyapp.md) | Developers & Archivists | Historical specifications and architectural lineage for legacy standalone prototype scripts (Version 1 & Version 2). | [Legacy Prototypes](https://byandell.github.io/MafADiscovery/docs/shinyapp.html) |

---

## 2. Core Repository Reference Documents

Outside of this `guides/` directory, the following primary documentation files are maintained at key repository locations:

* [**`../README.md`**](../README.md): Primary repository landing page, quick start instructions, environment setup, and feature highlights.
* [**`../DEVELOPER.md`**](../DEVELOPER.md): Master technical developer manual detailing data layer schemas, GRCm39 coordinates, reactive state architecture, and coding standards.
* [**`../SHINY_APP/README.md`**](../SHINY_APP/README.md): Dedicated user and collaborator manual for the production Shiny application in `SHINY_APP/`.
* [**`../AGENTS.md`**](../AGENTS.md): Coding guidelines and governance rules for AI assistants working in this repository.

---

## 3. In-App Modal Integration

The application in [`../SHINY_APP/app.R`](../SHINY_APP/app.R) loads guides directly into interactive modal dialogs via `render_markdown_file()`:

* **ℹ️ Help Modal (`input$show_help`)**: Reads [`user_guide.md`](user_guide.md) to orient researchers on plot axes, symbols, and navigation.
* **📖 Dev Guide Modal (`input$show_dev_guide`)**: Reads [`developer_guide.md`](developer_guide.md) providing direct hyperlinks to web documentation and source code.

Path resolution searches local and parent paths so documentation resolves cleanly during local R sessions, containerized deployments, and Shinylive WebAssembly builds.

---

## 4. Documentation Rendering Pipeline

All Markdown files in this directory are compiled into standalone, styled HTML pages using [`../render_docs.R`](../render_docs.R):

```bash
Rscript render_docs.R
```

The resulting pages are placed into `docs/` and deployed to GitHub Pages via `.github/workflows/deploy-shinylive.yaml`, serving live documentation at [https://byandell.github.io/MafADiscovery/docs/](https://byandell.github.io/MafADiscovery/docs/).
