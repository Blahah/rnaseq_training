### Quality control of the expression data

RNA-seq is a complex set of protocols, and there are many potential sources of bias and error that can affect the results. Whenever possible, we should try to detect these.

In this step, we will perform some basic quality control analysis of the expression quantification results.

An R script is provided at `./analysis/02_quality_control.R`. It:

- loads the merged data tables we created in the last step
- produces a plot that enables comparing the expression distribution between samples
- produces a plot that helps inspect the within- and between-condition correlation of expression estimates
- saves both plots

Open the script `./analysis/02_quality_control.R` in your RStudio. The easiest way to do this is to look at the `Files` tab in the bottom-right pane, and then click the `analysis` directory, followed by the `02_quality_control.R` script.

**Read through the script comments and code**. Execute the code one piece at a time, making sure each piece runs successfully and thinking about what it does.

**Look carefully at each plot**. Think about what it tells you. What would you see in the ideal scenario? What might indicate a problem?

When you've finished go to the terminal and run:

```bash
ls ./plots/quality_control
```

You should see two new files: `stacked_TPM_histograms.png` and `between_sample_correlation_heatmap.png`.

Now we'll move on to differential expression.
