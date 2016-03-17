# In this final script we will run the differential expression analysis
# We go from expression estimates to statistical inference about whether
# each transcript is differentially expressed, and then filter the
# transcripts to see only the highest-confidence DE calls.

# load the merged count data, rounding counts to integers as we go
est_counts <- round(read.table('data/merged/est_counts.tsv'))

# load the DESeq2 package we will use for differential expression testing
source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")

# create a data frame describing the condition for each experiment
split_names <- strsplit(colnames(est_counts), '.rep.')
conditions <- data.frame(
  condition = unlist(lapply(split_names, function(x) { x[1] }))
)
rownames(conditions) <- names(est_counts)

# prepare the DESeqDataSet from our experimental design and data
dds <- DESeqDataSetFromMatrix(countData = est_counts,
                              colData = conditions,
                              design = ~ condition)

# run DESeq2 on the dataset
dds <- DESeq(dds)

# extract the results
res <- results(dds)

# take a look!
print(res)

# summarize some basic tallies
summary(res)

# how man adjusted p-values were less than 0.1?
sum(res$padj < 0.1, na.rm=TRUE)

# the results function has a lot of features that can help
# look in more detail at the results
# for example, you could filter using a 0.05 alpha rather than 0.1
res05 <- results(dds, alpha=0.05)
summary(res05)

# you can use the plotCounts function to extract data for a specific transcript
# and then plot that data using ggplot2
d <- plotCounts(dds, gene=which.min(res$padj), intgroup="condition",
                returnData=TRUE)
library("ggplot2")
ggplot(d, aes(x=condition, y=count)) +
  geom_point(position=position_jitter(w=0.1,h=0)) +
  scale_y_log10(breaks=c(25,100,400))
