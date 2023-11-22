library(seqinr)
library(tidyverse)

# This script reads in the FASTA file for the T6SS genes
# and creates a BED format file listing each gene in the database.

# Import T6SS sequences.
seqsraw <- read.fasta("T6SS_CTD.fasta")
seqs <- sapply(seqsraw, function(x) getSequence(x, as.string=TRUE))
# Convert to a data frame.
dataSeqs <- data.frame(names(seqs), unname(unlist(seqs)))
colnames(dataSeqs) <- c("sample","seq")

# Build a BED file with one element per FASTA sequence.
dataSeqs <- dataSeqs %>%
  mutate(seqLength=nchar(seq))

# Generate a BED-format dataframe.
# The fields are:
# 1. chrom - chromosome name (in this case, the sequence name)
# 2. chromStart - starting position of the feature (in this case, 0, the first base)
# 3. chromEnd - ending position of the feature (in this case, the gene length)
dataBED <- dataSeqs %>%
  mutate(chromStart=0) %>%
  dplyr::select(sample, chromStart, seqLength) %>%
  dplyr::rename(chrom=sample, chromEnd=seqLength)

# Export the BED annotations as a table.
write.table(dataBED, "T6SS_CTD.bed",
            col.names=FALSE, quote=FALSE, row.names=FALSE, eol="\n", sep="\t")
