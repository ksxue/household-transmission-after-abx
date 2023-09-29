# This script loads the fixed differences dataframe and parses out information
# about households and timepoints for downstream analyses.

library(tidyverse)
library(data.table)

# Import the fixed differences dataframe.
dataFixedDiffs <- fread("workflow/analysis/identifyContaminatedSamples/out/fixedDiffs.txt.gz",
                        header=TRUE, data.table=FALSE)
# Parse the subjects and households from the information about fixed differences.
dataFixedDiffs <- dataFixedDiffs %>%
  mutate(subject1=substr(sample1,1,3), subject2=substr(sample2,1,3),
         hh1=substr(sample1,1,2), hh2=substr(sample2,1,2),
         timepoint1=as.integer(substr(sample1,5,7)),
         timepoint2=as.integer(substr(sample2,5,7)),
         type=ifelse(subject1==subject2,"sameSubject",
                     ifelse(hh1==hh2,"sameHousehold","diffHousehold")))

# Import the fixed differences summary.
dataFixedDiffsSummary <- fread("workflow/analysis/identifyContaminatedSamples/out/fixedDiffsSummary.txt",
                               header=TRUE, data.table=FALSE)

# Identify the species for which there is at least one high-coverage sample
# from at least two subjects in the same household.
highCovSpecies <- (dataFixedDiffsSummary %>%
                     filter(numSubjects>2, numWithinHouseholdPairs>1))$species
