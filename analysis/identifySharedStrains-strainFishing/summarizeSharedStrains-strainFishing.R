# This script takes the list of sharing annotations for within-household samples
# based on strain fishing, and it summarizes these sharing annotations to output
# a list of species that are shared and not shared for each household at key timepoints.

library(tidyverse)

# Import information on study metadata.
source("workflow/analysis/background/background.R")
# Set an output directory.
OUTDIR <- "workflow/analysis/identifySharedStrains-strainFishing/out/"

# Import the list of strain sharing annotations.
dataStrainFishingShared <- read.table(paste0(OUTDIR, "strainFishing-sharingAnnotated.txt.gz"),
                                      header=TRUE, stringsAsFactors = FALSE)

# Extract the comparison at the study initial timepoint, which is the baseline.
# This summary includes all subjects in the study, not only the ones in study arm X.
dataInitial <- dataStrainFishingShared %>%
  filter(bait %in% samplesInitial, sample %in% samplesInitial,
         type=="sameHousehold") %>%
  dplyr::select(species, bait, sample, baitSubject, sampleSubject, shared)
# Combine the strain sharing calls from the two directionalities of strain fishing.
# If the two "directions" of calls disagree, indicate that a strain is shared,
# since it is possible for a majority strain in one subject
# to be a minority strain in another, so that strain sharing is not bidirectional.
dataInitialSharing <- dataInitial %>%
  mutate(subject1=ifelse(baitSubject<sampleSubject, baitSubject, sampleSubject),
         subject2=ifelse(baitSubject<sampleSubject, sampleSubject, baitSubject)) %>%
  group_by(species, subject1, subject2) %>%
  summarize(shared=sum(shared)/n()) %>%
  mutate(sharedStrainFishing=!(shared==0))
# Export the initial strain sharing calls.
write.table(dataInitialSharing %>% dplyr::select(species, subject1, subject2, sharedStrainFishing), 
            paste0(OUTDIR, "sharedStrains-initial-all.txt"),
            row.names=FALSE, quote=FALSE)

# Extract the comparison across various key timepoints in the study.
# Include only subjects in study arm X for this analysis.
dataKeyTimepoints <- dataStrainFishingShared %>%
  filter((bait %in% samplesInitial & sample %in% samplesInitial) |
         (bait %in% samplesPreAbx & sample %in% samplesPreAbx) |
         (bait %in% samplesPostAbx & sample %in% samplesPostAbx) |
         (bait %in% samplesFinal & sample %in% samplesFinal)) %>%
  mutate(keyTimepoint=ifelse(baitTimepoint<29, "initial",
                             ifelse(baitTimepoint<34, "pre-abx",
                                    ifelse(baitTimepoint<60, "post-abx", "final")))) %>%
  dplyr::select(species, bait, sample, baitSubject, sampleSubject, baitHh, sampleHh, keyTimepoint, shared)
# Combine the strain sharing calls from the two directionalities of strain fishing.
# If the two "directions" of calls disagree, indicate that a strain is shared,
# since it is possible for a majority strain in one subject
# to be a minority strain in another, so that strain sharing is not bidirectional.
dataKeyTimepointsSharing <- dataKeyTimepoints %>%
  mutate(sample1=ifelse(bait<sample, bait, sample),
         sample2=ifelse(bait<sample, sample, bait),
         subject1=substr(sample1,1,3), subject2=substr(sample2,1,3),
         hh1=substr(sample1,1,2), hh2=substr(sample2,1,2)) %>%
  group_by(species, sample1, sample2, subject1, subject2, hh1, hh2, keyTimepoint) %>%
  summarize(shared=sum(shared)/n()) %>%
  mutate(sharedStrainFishing=!(shared==0))
# Export the initial strain sharing calls.
write.table(dataKeyTimepointsSharing %>% 
              dplyr::select(species, sample1, sample2, subject1, subject2, hh1, hh2, keyTimepoint, sharedStrainFishing), 
            paste0(OUTDIR, "sharedStrains-keyTimepoints-X.txt"),
            row.names=FALSE, quote=FALSE)


# Summarize the distribution of strain sharing in each cohabiting pair at each timepoint.
# Note that this only includes subjects from study arm X.
dataKeyTimepointsHouseholdSummary <- dataKeyTimepointsSharing %>%
  group_by(keyTimepoint, subject1, subject2) %>% dplyr::count(shared) %>%
  mutate(pctPopulations=n/sum(n))
write.table(dataKeyTimepointsHouseholdSummary, 
            paste0(OUTDIR, "sharedStrainsSummary-byHousehold-keyTimepoints-X.txt"),
            row.names=FALSE, quote=FALSE)
# Summarize the number of strains that are shared and not shared for each species.
dataKeyTimepointsSpeciesSummary <- dataKeyTimepointsSharing %>%
  group_by(species, keyTimepoint) %>% dplyr::count(shared) %>%
  mutate(pctPopulations=n/sum(n))
write.table(dataKeyTimepointsSpeciesSummary, 
            paste0(OUTDIR, "sharedStrainsSummary-bySpecies-keyTimepoints-X.txt"),
            row.names=FALSE, quote=FALSE)

