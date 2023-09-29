# This script consolidates the strain sharing calls
# from the fixed differences and strain fishing methods for each sample pair.
# For the strain fishing method, it also consolidates both "directions" of calls.
# It outputs a single table that includes each species-samplePair combination
# and lists whether strain sharing is identified for fixed differences or strain fishing.
# This script focuses only on calls from pairs of samples from the initial timepoint
# in study arm X. This includes both cohabiting and non-cohabiting pairs.
library(tidyverse)

# Import background metadata on the study.
source("workflow/analysis/background/background.R")
# Set the output directory.
OUTDIR <- "workflow/analysis/compareSharedStrainCalls/out/"

# Import the list of shared strain calls based on fixed differences.
dataFixedDiffsSharing <- 
  read.table("workflow/analysis/identifySharedStrains-fixedDiffs/out/fixedDiffs-initialX.txt.gz",
             header=TRUE, stringsAsFactors = FALSE)
# Trim down the number of relevant columns.
dataFixedDiffsSharing <- dataFixedDiffsSharing %>%
  dplyr::select(species, sample1, sample2, 
                subject1, subject2, hh1, hh2, timepoint1, timepoint2, shared) %>%
  dplyr::rename(sharedFixedDiffs=shared)

# Import the list of shared strain calls based on strain fishing.
dataStrainFishingSharing <-
  read.table("workflow/analysis/identifySharedStrains-strainFishing/out/strainFishing-sharingAnnotated-initialX.txt.gz",
             header=TRUE, stringsAsFactors = FALSE)
# Consolidate the two "directions" of shared strain calls from strain fishing.
# If the two "directions" of calls disagree, indicate that a strain is shared,
# since it is possible for a majority strain in one subject
# to be a minority strain in another, so that strain sharing is not bidirectional.
dataStrainFishingSharingConsolidated <- dataStrainFishingSharing %>%
  mutate(sample1=ifelse(bait<sample, bait, sample),
         sample2=ifelse(bait<sample, sample, bait),
         subject1=substr(sample1,1,3), subject2=substr(sample2,1,3),
         hh1=substr(sample1,1,2), hh2=substr(sample2,1,2),
         timepoint1=as.numeric(substr(sample1,5,7)), timepoint2=as.numeric(substr(sample2,5,7))) %>%
  group_by(species, sample1, sample2, subject1, subject2, hh1, hh2, timepoint1, timepoint2) %>%
  summarize(shared=sum(shared)/n()) %>%
  mutate(sharedStrainFishing=!(shared==0)) %>%
  dplyr::select(-shared)

# Combine the two dataframes of sample calls from the two methods.
dataSharedStrainsConsolidated <- 
  full_join(dataFixedDiffsSharing, dataStrainFishingSharingConsolidated,
            by=c("species","sample1","sample2","subject1","subject2",
                 "hh1","hh2","timepoint1","timepoint2")) %>%
  arrange(sample1, sample2, species)
# Annotate calls as concordant or discordant.
dataSharedStrainsConsolidated <- dataSharedStrainsConsolidated %>%
  mutate(concordant=(sharedFixedDiffs==sharedStrainFishing))

# Annotate consensus calls of strain sharing based on fixed differences and strain fishing.
# Use five methods of annotating this consensus:
# 1. conservative - all discordant cases are marked as false (no strain sharing)
# 2. permissive - all discordant cases are marked as true (strain sharing)
# 3. fixedDiffs - in discordant cases, use the call from fixed differences
# 4. strainFishing - in discordant cases, use the call from strain fishing
# 5. NA - in discordant cases, do not make a call
dataSharedStrainsConsolidated <- dataSharedStrainsConsolidated %>%
  mutate(consensusConservative=
           ifelse(is.na(concordant), 
                  ifelse(is.na(sharedFixedDiffs), sharedStrainFishing, sharedFixedDiffs),
                  ifelse(concordant, sharedStrainFishing, FALSE)),
         consensusPermissive=
           ifelse(is.na(concordant), 
                  ifelse(is.na(sharedFixedDiffs), sharedStrainFishing, sharedFixedDiffs),
                  ifelse(concordant, sharedStrainFishing, TRUE)),
         consensusFixedDiffs=
           ifelse(is.na(concordant), 
                  ifelse(is.na(sharedFixedDiffs), sharedStrainFishing, sharedFixedDiffs),
                  ifelse(concordant, sharedStrainFishing, sharedFixedDiffs)),
         consensusStrainFishing=
           ifelse(is.na(concordant), 
                  ifelse(is.na(sharedFixedDiffs), sharedStrainFishing, sharedFixedDiffs),
                  ifelse(concordant, sharedStrainFishing, sharedStrainFishing)),
         consensusNA=
           ifelse(is.na(concordant), 
                  ifelse(is.na(sharedFixedDiffs), sharedStrainFishing, sharedFixedDiffs),
                  ifelse(concordant, sharedStrainFishing, NA)))

# Annotate samples that derive from key timepoints.
dataSharedStrainsConsolidated <- dataSharedStrainsConsolidated %>%
  mutate(keyTimepoint1=ifelse(sample1 %in% samplesInitial, "initial",
                       ifelse(sample1 %in% samplesPreAbx, "pre-abx",
                       ifelse(sample1 %in% samplesPostAbx, "post-abx",
                       ifelse(sample1 %in% samplesFinal, "final", NA)))),
         keyTimepoint2=ifelse(sample2 %in% samplesInitial, "initial",
                       ifelse(sample2 %in% samplesPreAbx, "pre-abx",
                       ifelse(sample2 %in% samplesPostAbx, "post-abx",
                       ifelse(sample2 %in% samplesFinal, "final", NA)))))

# Export the consolidated set of strain sharing calls.
write.table(dataSharedStrainsConsolidated,
            paste0(OUTDIR, "strainSharingCalls-consolidated-initialX.txt"),
            row.names=FALSE, quote=FALSE)
