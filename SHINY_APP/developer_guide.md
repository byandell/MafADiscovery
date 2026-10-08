# MafA Discovery: Developer Guide & Architecture

This explorer is documented across four core architectural modules maintained in the repository:

### Core Documentation Modules

1. **[Developer Guide (DEVELOPER.md)](https://github.com/byandell/MafADiscovery/blob/main/DEVELOPER.md)** ([Site Link](DEVELOPER.md))  
   Master technical guide covering repository architecture, mouse GRCm39 coordinate mapping, reactive lifecycle, coding standards, and deployment rules.

2. **[Legacy Prototypes (shinyapp.md)](https://github.com/byandell/MafADiscovery/blob/main/shinyapp.md)** ([Site Link](shinyapp.md))  
   Technical specifications and lineage for Version 1 (`MafA_Discovery_App.R`) and Version 2 (`MafA_Discovery_App_v2.R`) standalone prototype scripts.

3. **[Publishing & Deployment (publishapp.md)](https://github.com/byandell/MafADiscovery/blob/main/publishapp.md)** ([Site Link](publishapp.md))  
   Shinylive (webR) static WebAssembly deployment architecture, performance benchmarks, and automated GitHub Actions CI configuration.

4. **[App Redesign & Reactivity (redesign.md)](https://github.com/byandell/MafADiscovery/blob/main/redesign.md)** ([Site Link](redesign.md))  
   Dynamic view mode lifecycle (`Genome-Wide`, `Chromosome`, `QTL Region`, `Locus Zoom`), discrete filter ergonomics, coordinate auto-scaling, and total state reset.

5. **[F2 Glycemic QTL Integration (qtlanalysis.md)](https://github.com/byandell/MafADiscovery/blob/main/qtlanalysis.md)** ([Site Link](qtlanalysis.md))  
   F2 glycemic QTL integration, automated 95% CI bounding (`ci.low` to `ci.high`), visual interval highlights, and interactive QTL reference table.

---

### Repository Links

* **GitHub Repository:** [byandell/MafADiscovery](https://github.com/byandell/MafADiscovery)
* **Application Source:** [`SHINY_APP/app.R`](https://github.com/byandell/MafADiscovery/blob/main/SHINY_APP/app.R)
