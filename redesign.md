# Redesign Ideas

## Ideas to improve reactivity and navigation

* Hide "Select Chromosome:" when in "View Mode: Genome-Wide"
* "View Mode: Locus Zoom" 
  * If "Search Gene Symbol:" is empty, main panel is blank.
  * Don't show "Locus Zoom" option until Gene Symbol is non-empty.
  * Once a "Gene Symbol" is selected, it appears.
  * If you then blank out "Search Gene Symbol", the last searched gene symbol is still selected.
* "Locus Window (Kbp):" this is only relevant when in "Locus Zoom"
  * Might be useful to have only selected set of values, say 20,50,100,200,500,1000,2000
* When "Select Gene Symbol"
  * Top panel initially has whole genome but then goes to its chromosome;
chromosome needs to be named (only says "Chromosome").
* "Visible Categories:" is redundant given plotly ability to click on labels.
* "Min Peak Score:" consider selected values as for "Locus Window".

## Challenges with Interpreting

* Need explanation of both panels on "Locus Zoom".
This could be pop-up or pull-down.
* Reactive hiding and showing of side panel would simplify presentation
(see ideas above).

