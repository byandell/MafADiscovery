# Mark Notes

**Prompt:**
Develop implementation plan for features and modifications to the app described below.

## New Features and Modifications

New feature: integrate loci identified in the F2 study, as listed in
[Top_glycemic_QTL_for_sex_additive_analysis.csv](SHINY_APPs/Top_glycemic_QTL_for_sex_additive_analysis.csv).

We want to show this table and connect these QTL regions to other features.
For example, to zoom into the QTL on Chr16, I select view mode=chromosome, select chromosome=11, and then use the zoom tool for ~50 to 100Mbp. I then click on one of the triangles to highlight, yielding this view, where three genes are DE, all going up in the SJL backcrossed mice and proximal to a MAFA peak with SNPs.

Jumping from this locus on Chr16 to our QTL on 13, requires select chromosome 13 and then zoom for ~0 to 50Mbp, yielding this after clicking on a orange diamond (coding variant).

Since this type of query will be the focus on the app when integrating the QTL from the F2s, what would it take to include the QTL as another option for view mode? We could offer Chr and the corresponding CI for low and high as boundaries that will be displayed.

Finally, deselect “Shared” from the visible categories as the default.
