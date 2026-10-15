# About MafA Discovery

An interactive genomic explorer for MafA transcription factor binding peaks, strain-divergent variants, differential gene expression, and F2 glycemic QTL loci in C57BL/6J and SJL/J mouse models.

* [Research Overview](#research-overview)
* [Methods](#methods)
* [Integrated Data Layers](#integrated-data-layers)
* [Collaborative Consortium](#collaborative-consortium)
* [Software & Implementation](#software--implementation)

---

## Research Overview

MafA is a critical basic leucine zipper (bZIP) transcription factor that regulates pancreatic $\beta$-cell maturation, insulin gene transcription, and glucose homeostasis. Genetic variation between inbred mouse strains alters transcription factor binding affinity, chromatin accessibility, and downstream transcriptional networks, predisposing certain strains to $\beta$-cell dysfunction and diabetes.

The **MafA Discovery Integrated Genomic Explorer** was developed to bridge molecular binding events with organismal phenotypic outcomes:

* Linking MafA ChIP-seq / CUT&RUN binding peaks to proximal strain-divergent genetic variants between **C57BL/6J** (diabetes-resistant baseline) and **SJL/J** (diabetes-susceptible).
* Connecting MafA chromatin occupancy to downstream differential gene expression (DEGs) in islet tissue.
* Mapping candidate binding peaks and high-impact coding variants directly into 95% confidence intervals of **F2 glycemic QTL loci**.

---

## Methods

Genetic analysis was performed using quantitative trait loci (QTL) methods implmented in R/qtl (Broman et al., 2003, 2019). This involved genome scans of GTT measurements,
including summary traits based on AUC (area under the curve) and further
meta-trates over the multiple weeks of assays. Peaks and 95% approximate
confidence intervals were reported.

The peak scores were developed ...

### Methods References

* [R/qtl2 QTL analysis for high-dimensional data and complex crosses (Karl Broman)](https://kbroman.org/qtl2/)
* [Broman K. W., Wu H., Sen S., Churchill G. A., 2003. R/qtl: QTL mapping in experimental crosses. Bioinformatics 19: 889–890.](https://doi.org/10.1093/bioinformatics/btg112)
* [Broman KW, Gatti DM, Simecek P, Furlotte NA, Prins P, Sen Ś, Yandell BS, Churchill GA (2019) R/qtl2: software for mapping quantitative trait loci with high-dimensional data and multi-parent populations. Genetics 211:495-502](https://doi.org/10.1534/genetics.118.301595)

## Integrated Data Layers

All genomic coordinates across datasets are aligned to the mouse **GRCm39 (mm39)** reference assembly:

1. **MafA Binding Peaks (`MafA_Peaks_with_SNPs_v3.csv`)**:
   * Peak intervals, genomic midpoints, and ChIP/CUT&RUN binding strength scores.
   * Total strain-divergent variant counts within each binding peak interval (`variant_count`).
2. **Differential Gene Expression (`Master_DEG_Strain_Comparison_v3.csv`)**:
   * RNA-seq differential expression between C57BL/6J and SJL/J islets under baseline and perturbational conditions.
   * Categorized into `Shared`, `C57_Specific`, `SJL_Specific`, and `Discordant` responses with log2 fold-changes.
3. **Mouse Genomic Backbone (`mouse_genes_mm39_v3.csv`)**:
   * Comprehensive Ensembl GRCm39 gene models with coordinates, strand orientations, and gene biotypes.
4. **Prioritized Coding Polymorphisms (`B6_SJL_prioritized_protein_coding_SNPs.csv`)**:
   * Protein-coding sequence variants between B6 and SJL classified by impact severity (`HIGH` and `MODERATE`).
   * Annotated with evolutionary conservation scores (`phastCons`), consequence types (`csq`), and amino acid substitutions (`aa_change`).
5. **F2 Glycemic QTL Loci & Scans (`QTLresults/Top_glycemic_QTL_for_sex_additive_analysis.csv`)**:
   * Top quantitative trait loci from the B6 $\times$ SJL F2 intercross study across glycemic traits (e.g., AUC, glucose response).
   * Annotated with peak marker coordinates, LOD scores, 95% Bayesian confidence intervals (`ci.low` to `ci.high`), additive effect models (`BB`, `BS`, `SS`), and corresponding trait scan plots.

---

## Collaborative Consortium

This investigation is a joint collaboration between:

* **Vanderbilt University Medical Center**:
  * Department of Molecular Physiology and Biophysics
  * Jee Yeon Cha (PI) and Mallory Maurer
* **University of Wisconsin-Madison**:
  * Department of Statistics — Brian Yandell
  * Department of Biochemistry — Mark Keller

---

## Software & Implementation

* **Production Application**: Built with R Shiny, `data.table`, `ggplot2`, and `plotly`.
* **Serverless WebAssembly (WASM)**: Deployed serverlessly via [Shinylive for R (webR)](https://github.com/posit-dev/shinylive-r) on GitHub Pages, requiring zero server maintenance.
* **Open Source Repository**: Source code and documentation are available on GitHub at [byandell/MafADiscovery](https://github.com/byandell/MafADiscovery).
