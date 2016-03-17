# In this script we take the expression quantification tables
# output by salmon, one per sample, and we load them.
#
# We collect the measures we are interested in into one data frame
# per measure, with one row per transcript and one column per sample.
#
# Finally we use the metadata for the experiment to name the columns
# in out collected data frames, and then we write out the data frames
# to a new 'merged' directory.

# we want a separate data frame for each of:
# - expression (TPM)
# - expression (estimated counts)
# - effective length
# we create the variables but leave them NULL for now
# they'll get filled in as we load the data
tpm <- NULL
est_counts <- NULL
eff_length <- NULL

# load the files for each sample in a loop
files <- Sys.glob('data/salmon/quant_SRR313673*/quant.sf')
for (file in files) {
  df <- read.table(file, header=TRUE, stringsAsFactors=FALSE)

  first <- is.null(tpm)

  # store the different columns
  tpm <- cbind(tpm, df[, 'TPM'])
  est_counts <- cbind(est_counts, df[, 'NumReads'])
  eff_length <- cbind(eff_length, df[, 'EffectiveLength'])

  # make sure the collected results have the trancript IDs
  if (first) {
    rownames(tpm) <- df$Name
    rownames(est_counts) <- df$Name
    rownames(eff_length) <- df$Name
  }

  rm(df)
}

# match up the filenames with the sample and replicate information
# load the sample information
samples <- read.table('data/metadata/E-GEOD-77428.sdrf.txt',
                      header=TRUE, sep="\t")

# extract the SRR identifier from each result file path
srrs <- gsub('.+/quant_(SRR[0-9]+)/.*', '\\1', files)

# extract the sample names from the samples table in the same order
# as our files were loaded (so, the same order as the columns in our results)
sample_names <- as.character(samples[match(srrs, samples$Scan.Name), 'Comment..Sample_title.'])

# now we can name the columns in our results!
colnames(tpm) <- sample_names
colnames(est_counts) <- sample_names
colnames(eff_length) <- sample_names

# and finally, we save the merged data
write.table(tpm, 'data/merged/tpm.tsv', sep="\t")
write.table(est_counts, 'data/merged/est_counts.tsv', sep="\t")
write.table(eff_length, 'data/merged/eff_length.tsv', sep="\t")

# now head back to the browser to read instructions for the next stage of analysis
