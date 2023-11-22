# This script consolidates all samples that should be excluded from
# downstream analyses.
# These samples may have evidence of cross-contamination,
# low sequencing coverage, or suspicion of sample swaps.

library(tidyverse)
library(data.table)

# Import the background list of samples.
source("workflow/analysis/background/background.R")

# Import the list of samples with evidence of cross-contamination.
# First test to make sure that it exists.
contaminatedSamples <- c()
if(file.exists("workflow/analysis/identifyContaminatedSamples/out/contaminatedSamples.txt")){
  contaminatedSamples <- (read.table("workflow/analysis/identifyContaminatedSamples/out/contaminatedSamples.txt",
                                     header=TRUE, stringsAsFactors = FALSE))$sample
}


# Exclude all follow-up samples from household H from downstream analyses.
# The species abundance plots for these subjects indicate that
# there was some swapping of subject IDs during follow-up sampling.
# The subjects also indicated this at the time of sample pickups.
swappedSamples <- (samplesRaw %>%
  filter(hh=="XH", timepoint>75))$sample

# Aggregate all samples that should be excluded from downstream analyses.
# Also exclude sample WWW-000, which came from an unknown subject.
sampleBlacklist <- as.data.frame(sort(unique(c(contaminatedSamples, swappedSamples, "WWW-000"))))
colnames(sampleBlacklist) <- c("sample")
# Export list of blacklisted samples.
write.table(sampleBlacklist, "workflow/analysis/background/out/sampleBlacklist.txt",
            quote=FALSE, row.names=FALSE)
