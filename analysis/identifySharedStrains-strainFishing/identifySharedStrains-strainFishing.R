# This script takes in the filtered set of strain fishing calculations.
# It calculates the percent of private sites detected in each bait-sample pair
# and uses this percent to make inferences about which strains are shared
# by cohabiting subjects in the study.

library(tidyverse)

# Set the output directory for these analyses.
outdir <- "workflow/analysis/identifySharedStrains-strainFishing/out/"

# Load the filtered strain fishing dataframe.
source("workflow/analysis/identifySharedStrains-strainFishing/loadStrainFishingFiltered.R")
# Load the background metadata.
source("workflow/analysis/background/background.R")

# Annotate bait-sample pairs as shared or not shared based on the percent of available
# private SNPs detected in the sample.
dataStrainFishing <- dataStrainFishing %>%
  mutate(shared=(pctSitesDetected>0.75))

# Export the number of private SNPs per bait sample for each species.
numPrivateSNPs <- dataStrainFishing %>%
  filter(bait==sample) %>%
  dplyr::select(species, bait, numSitesAvailable) %>%
  dplyr::rename(numPrivateSNPsTotal=numSitesAvailable)
write.table(numPrivateSNPs, 
            paste0(outdir,"numPrivateSNPsPerBait.txt"),
            quote=FALSE, row.names=FALSE)

# Export the full set of sharing annotations for samples from the same subject.
strainSharingAnnotationsSameSubject <- dataStrainFishing %>%
  filter(type=="sameSubject", bait!=sample)
write.table(strainSharingAnnotationsSameSubject, 
            gzfile(paste0(outdir,"strainFishing-sameSubject.txt.gz")),
            quote=FALSE, row.names=FALSE)

# Export the full set of sharing annotations for samples from the same household.
strainSharingAnnotations <- dataStrainFishing %>%
  filter(type=="sameHousehold")
write.table(strainSharingAnnotations, 
            gzfile(paste0(outdir,"strainFishing-sharingAnnotated.txt.gz")),
            quote=FALSE, row.names=FALSE)

# Export the sharing annotations for all initial samples from study arm X.
strainSharingAnnotationsInitialX <- dataStrainFishing %>% 
  filter(bait %in% samplesInitial, sample %in% samplesInitial,
         baitSubject %in% subjectsX, sampleSubject %in% subjectsX,
         bait!=sample)
write.table(strainSharingAnnotationsInitialX, 
            gzfile(paste0(outdir,"strainFishing-sharingAnnotated-initialX.txt.gz")),
            quote=FALSE, row.names=FALSE)

# Export the full SNP frequency and strain fishing information for the bait samples,
# which provides information about the initial strains identified in each population.
baitStrainFrequencies <- dataStrainFishing %>%
  filter(bait==sample)
write.table(baitStrainFrequencies, paste0(outdir, "baitStrains.txt"),
            quote=FALSE, row.names=FALSE)

# Export the list of "cousin" strains; that is, populations from different households
# that have high similarity.
cousins <- dataStrainFishing %>%
  filter(type=="diffHousehold", shared)
write.table(cousins, 
            gzfile(paste0(outdir,"strainFishing-cousins.txt.gz")),
            quote=FALSE, row.names=FALSE)
# Calculate the percentage of between-household comparisons that exceed the
# 75% threshold, i.e. the percentage of cousin comparisons.
pctCousins <- dataStrainFishing %>%
  group_by(species) %>%
  summarize(numCousins=sum(type=="diffHousehold" & shared),
            numBtwnHousehold=sum(type=="diffHousehold"),
            pctCousins=numCousins/numBtwnHousehold)
write.table(pctCousins, 
            gzfile(paste0(outdir,"strainFishing-pctCousins.txt")),
            quote=FALSE, row.names=FALSE)

# Export all comparisons for a small set of sample species.
write.table(dataStrainFishing %>%
              filter(species %in% c("Bacteroides_fragilis_54507",
                                    "Bifidobacterium_longum_57796",
                                    "Eubacterium_eligens_61678")),
            gzfile(paste0(outdir, "strainFishing-sampleSpecies.txt.gz")),
            quote=FALSE, row.names=FALSE)
