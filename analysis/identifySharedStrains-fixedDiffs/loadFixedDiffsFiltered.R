# This script loads the fixed differences dataframe and parses out information
# about households and timepoints for downstream analyses.
# It also reads in the list of samples with signatures of cross-contamination
# and removes them from the dataset.

# Import the fixed differences dataframe.
source("workflow/analysis/identifyContaminatedSamples/loadFixedDiffs.R")
# Import the list of blacklisted samples.
sampleBlacklist <- (read.table("workflow/analysis/background/out/sampleBlacklist.txt",
                              header=TRUE, stringsAsFactors = FALSE))$sample

# Filter the fixed differences dataframe,
# removing all pairs that contain samples with evidence of cross-contamination.
dataFixedDiffs <- dataFixedDiffs %>%
  filter(!(sample1 %in% sampleBlacklist),
         !(sample2 %in% sampleBlacklist))
# Retain only pairs from high-coverage species
# for which there is at least one high-coverage sample
# from at least two subjects in the same household.
dataFixedDiffs <- dataFixedDiffs %>%
  filter(species %in% highCovSpecies)
