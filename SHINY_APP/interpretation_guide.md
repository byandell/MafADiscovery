# How to Interpret Panels & Navigate

### 1. Top Panel: Manhattan Macro View
Displays MafA binding peaks, prioritized coding SNPs, and F2 glycemic QTL intervals:

* **Primary Y-Axis (Left):** MafA Peak Score reflecting ChIP/CUT&RUN binding strength.
* **Point Shapes:** Differentially expressed gene (DEG) direction (▲ UP, ▼ DOWN in diabetes/perturbation).
* **Point Colors:** Strain-specificity category:
  * <span style="color:#228B22; font-weight:bold;">Shared</span> (green)
  * <span style="color:#0000CC; font-weight:bold;">C57_Specific</span> (blue)
  * <span style="color:#CC0000; font-weight:bold;">SJL_Specific</span> (red)
  * <span style="color:#FF8C00; font-weight:bold;">Discordant</span> (orange)
* **Secondary Y-Axis & Diamonds:** Coding SNPs between C57BL/6J and SJL/J strains, scaled by evolutionary conservation score (`phastCons` 0–1). Pink/purple = HIGH impact, Gold = MODERATE impact.
* **QTL Interval Highlight:** Translucent blue shading and vertical dashed line show the F2 glycemic QTL confidence interval (95% CI) and peak marker position.
* **Interactivity:** Click any peak or SNP point to inspect its locus. Use "Visible Categories" to isolate strain-divergent peaks.

---

### 2. Bottom Panel: Locus Micro Schematic
Renders a high-resolution window (± Locus Window) centered on the active MafA binding peak:

* **MafA Peak:** Marked with a gold dashed line and highlighted region.
* **Gene Models:** Horizontal bars depict gene bodies. Arrows denote transcription start site (TSS) and orientation.
* **Coding SNPs:** Diamonds show exact positions of coding SNPs within exons. Outlined diamonds indicate `phastCons` ≥ 0.7.
* **Peak Deselection:** Click **"✕ Deselect Current Peak"** (located in the Manhattan Controls side panel immediately below the Hide Legend checkbox) to clear the active peak and close the locus view.

---

### 3. Navigation Modes
* **Genome-Wide:** Global linear coordinate view across all 21 mouse chromosomes.
* **Chromosome:** Focused view of an individual chromosome with absolute Mbp tick marks.
* **QTL Region:** Access the interactive F2 Glycemic QTL Overview Table or zoom into the 95% confidence interval boundaries (`ci.low` to `ci.high`) of a selected QTL, displaying its additive sex QTL LOD scan directly above the Manhattan plot.
* **Locus Zoom:** Fine-scale schematic centered on an active MafA binding peak.
