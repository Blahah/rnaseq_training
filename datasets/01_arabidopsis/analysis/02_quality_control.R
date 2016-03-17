# In this script we load the merged expression data
# and perform some quality control analysis

# load the merged TPM data
tpm <- read.table('data/merged/tpm.tsv')

# check that the data has loaded OK
# there should be 6 named columns and 40,223 rows
head(tpm)
dim(tpm)

# we need a specific version of ggplot2, so install the latest one now:
install.packages('ggplot2')

# create a directory in which we can save our quality control plots
# we turn off warnings because we expect the directory might already exist
# and don't want to see a warning if it does
dir.create('plots/quality_control', showWarnings = FALSE)

# generally, we expect every RNA sample in an experiment to have a similar
# expression distribution.
#
# between any two conditions, only a small number of transcripts should
# be expressed differently. And between replicates of a single condition
# we hope no transcripts will be expressed differently.
#
# let's visualise the distribution of each sample and stack
# them so we can compare easily
#
# first we need to convert the TPM data to long format so it will
# work nicely with ggplot2
library(reshape2)
tpm_long <- melt(tpm) # this will print a warning - it's fine to ignore it
names(tpm_long) <- c('sample', 'tpm')

# now we want to have separate columns describing the condition and replicate
condition_rep <- strsplit(as.character(tpm_long$sample), '.rep.', fixed=TRUE)
tpm_long$condition <- factor(unlist(lapply(condition_rep, function(x) { x[1] })))
tpm_long$replicate <- factor(unlist(lapply(condition_rep, function(x) { x[2] })))

# with the data in the right shape, we can now make our stacked
# distribution plot
library(ggplot2)
stacked_tpm_histograms <- ggplot(tpm_long, aes(x=tpm, fill=condition)) +
  geom_histogram(bins = 100) +
  scale_x_log10(limits = c(1e-04, 1e+05)) +
  facet_grid(condition+replicate~.) +
  labs(fill = "Genotype",
       y = 'Number of transcripts',
       x = 'Transcripts per Million (TPM) on log10 scale') +
  ggtitle('Expression distribution for each sample')

# display the plot
# this will show some warnings, but they are harmless
print(stacked_tpm_histograms)

# save the plot
# this will show some warnings, but they are harmless
ggsave(filename = 'plots/quality_control/stacked_TPM_histograms.png',
       plot = stacked_tpm_histograms, width = 10, height = 10)

# this plot helps us quickly check to see if any sample looks way off.
# If any distribution looks very different from the rest, that suggests
# something went wrong with that sample.
#
# If all the distributions look very different
# from one another, we may have poor consistency among the samples.
#
# If everything looks good, we might next look in more detail. We
# expect that the replicates of a sample should correlate very well
# with one another, and that the correlation within a sample should be
# better than the between sample correlation.
#
# Let's take a look at the pairwise correlation between all samples

# first we calculate the correlation coefficient between each pair of samples
# we use the Spearman rank correlation because it isn't affected too badly
# by outliers
tpm_cor <- melt(cor(tpm, method='spearman'))

# now we'll make a heatmap of the correlation data
between_sample_correlation_heatmap <-
  ggplot(data.frame(tpm_cor), aes(x = Var1, y = Var2, fill = value)) +
  geom_bin2d(stat = 'identity') +
  geom_text(aes(label = round(value, 3))) +
  labs(fill = 'Spearman\ncorrelation', x = 'Sample', y = 'Sample') +
  ggtitle('Pairwise Spearman correlation\nbetween all samples')

# display the plot
print(between_sample_correlation_heatmap)

# save the plot
ggsave(filename = 'plots/quality_control/between_sample_correlation_heatmap.png',
       plot = between_sample_correlation_heatmap, width = 10, height = 10)

# Here's what we're looking for:
# - the correlation between reps of a condition should be high (> 0.9 usually)
# - the correlation between conditions should also be high,
#   but less high than the within-condition correlation
# If those conditions are met - our data is looking good!
