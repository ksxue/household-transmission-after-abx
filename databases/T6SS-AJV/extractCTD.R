library(seqinr)
library(tidyverse)

# This script reads in the FASTA file for the GA2 effectors
# and extracts the CTD, corresponding to the last 400 bases of the sequence.

# Import GA2 effector sequences.
seqsraw <- read.fasta("GA2_effectors.fasta")
seqs <- sapply(seqsraw, function(x) getSequence(x, as.string=TRUE))
# Convert to a data frame.
dataSeqs <- data.frame(names(seqs), unname(unlist(seqs)))
colnames(dataSeqs) <- c("sample","seq")

# Extract the CTD for each sequence.
dataSeqs <- dataSeqs %>%
  mutate(seqLength=nchar(seq),
         CTD=toupper(substr(seq, seqLength-400+1, seqLength)),
         CTDLength=nchar(CTD),
         seqName=paste0(sample, "_CTD"))

# Export the CTDs as a FASTA file.
write.fasta(as.list(dataSeqs$CTD), dataSeqs$seqName,
            "GA2_effectors_CTD.fasta")


# This script reads in the FASTA file for the GA3 effectors
# and extracts the CTD, corresponding to the last 400 bases of the sequence.

# Import GA3 effector sequences.
seqsraw <- read.fasta("GA3_effectors.fasta")
seqs <- sapply(seqsraw, function(x) getSequence(x, as.string=TRUE))
# Convert to a data frame.
dataSeqs <- data.frame(names(seqs), unname(unlist(seqs)))
colnames(dataSeqs) <- c("sample","seq")

# Extract the CTD for each sequence.
dataSeqs <- dataSeqs %>%
  mutate(seqLength=nchar(seq),
         CTD=toupper(substr(seq, seqLength-400+1, seqLength)),
         CTDLength=nchar(CTD),
         seqName=paste0(sample, "_CTD"))

# Export the CTDs as a FASTA file.
write.fasta(as.list(dataSeqs$CTD), dataSeqs$seqName,
            "GA3_effectors_CTD.fasta")
