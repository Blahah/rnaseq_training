### Differential expression analysis

Finally, we come to the differential expression analysis step.

This step can be fairly simple - and in this example we will make it as simple as possible. If there is time, we can explore more detailed aspects of the analysis.


An R script is provided [here](http://rik.smith-unna.com/rnaseq_training/datasets/01_arabidopsis/analysis/03_differential_expression.R). Right-click the link and save it as `~/Desktop/RNA/analysis/03_differential_expression.R`

It:

- loads the merged data tab les we created in the last step
- produces a plot that enables comparing the expression distribution between samples
- produces a plot that helps inspect the within- and between-condition correlation of expression estimates
- saves both plots

Open the script `./analysis/03_differential_expression.R` in your RStudio. The easiest way to do this is to look at the `Files` tab in the bottom-right pane, and then click the `analysis` directory, followed by the `03_differential_expression.R` script.

**Read through the script comments and code**. Execute the code one piece at a time, making sure each piece runs successfully and thinking about what it does.

**Once you've reached the end of this script** you should congratulate yourself! You've just successfully completed a differential expression experiment!

However, this was just a taste of how to perform a full differential expression analysis. To explore what else you should think about, the following resources provide a great deal more detail:

- http://www-huber.embl.de/users/klaus/Teaching/DESeq2Predoc2014.html
- http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.pdf
