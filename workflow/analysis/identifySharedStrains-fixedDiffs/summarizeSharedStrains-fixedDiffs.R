# This script takes the annotated dataframe of fixed differences
# between pairs of samples of each high-coverage species.
# It generates several summaries of the number of strains shared at different
# key timepoints.

library(tidyverse)
library(data.table)

# Import background information on study design.
source("workflow/analysis/background/background.R")
# Set an output directory.
OUTDIR <- "workflow/analysis/identifySharedStrains-fixedDiffs/out/"

# Import the annotated summary of fixed differences.
# Note that this dataframe includes strains annotated as shared and not shared.
dataFixedDiffsShared <- read.table(paste0(OUTDIR,"fixedDiffs-sharingAnnotated.txt.gz"),
                                   header=TRUE, stringsAsFactors = FALSE)

# Extract the comparison at the study initial timepoint, which is the baseline.
# This summary includes all subjects in the study, not only the ones in study arm X.
dataInitial <- dataFixedDiffsShared %>%
  filter(sample1 %in% samplesInitial, sample2 %in% samplesInitial) %>%
  dplyr::select(sample1, sample2, subject1, subject2, hh1, hh2, species, shared)
write.table(dataInitial, paste0(OUTDIR, "sharedStrains-initial-all.txt"),
            row.names=FALSE, quote=FALSE)


# Extract the comparison across various key timepoints in the study.
# Include only subjects in study arm X for this analysis.
dataKeyTimepoints <- dataFixedDiffsShared %>%
  filter((sample1 %in% samplesInitial & sample2 %in% samplesInitial) |
         (sample1 %in% samplesPreAbx & sample2 %in% samplesPreAbx) |
         (sample1 %in% samplesPostAbx & sample2 %in% samplesPostAbx) |
         (sample1 %in% samplesFinal & sample2 %in% samplesFinal)) %>%
  filter(sample1 %in% samplesKeyTimepoints & sample2 %in% samplesKeyTimepoints) %>%
  mutate(keyTimepoint=ifelse(timepoint1<29, "initial",
                      ifelse(timepoint1<34, "pre-abx",
                      ifelse(timepoint1<60, "post-abx", "final")))) %>%
  dplyr::select(sample1, sample2, subject1, subject2, hh1, hh2, keyTimepoint, species, shared)
write.table(dataKeyTimepoints, paste0(OUTDIR, "sharedStrains-keyTimepoints-X.txt"),
            row.names=FALSE, quote=FALSE)

# Summarize the distribution of strain sharing in each cohabiting pair at each timepoint.
# Note that this only includes subjects from study arm X.
dataKeyTimepointsHouseholdSummary <- dataKeyTimepoints %>%
  group_by(keyTimepoint, subject1, subject2) %>% dplyr::count(shared) %>%
  mutate(pctPopulations=n/sum(n))
write.table(dataKeyTimepointsHouseholdSummary, 
            paste0(OUTDIR, "sharedStrainsSummary-byHousehold-keyTimepoints-X.txt"),
            row.names=FALSE, quote=FALSE)
# Summarize the number of strains that are shared and not shared for each species.
dataKeyTimepointsSpeciesSummary <- dataKeyTimepoints %>%
  group_by(species, keyTimepoint) %>% dplyr::count(shared) %>%
  mutate(pctPopulations=n/sum(n))
write.table(dataKeyTimepointsSpeciesSummary, 
            paste0(OUTDIR, "sharedStrainsSummary-bySpecies-keyTimepoints-X.txt"),
            row.names=FALSE, quote=FALSE)
