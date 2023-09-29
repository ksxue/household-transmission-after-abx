# This script imports the species abundances calculated by MIDAS
# and compares the species shared between every pair of subjects
# at the initial study timepoint in terms of their number and relative abundance.

library(tidyverse)
library(foreach)

# Import metadata on the samples and subjects.
source("workflow/analysis/background/background.R")
# Set the output directory.
OUTDIR <- "workflow/analysis/compareSpeciesAbundances/out/"
# Set the limit of detection for species detected through metagenomics.
LIMITOFDETECTION <- 0.001


# Import the list of species abundances.
source("workflow/analysis/generateSpeciesAbundances/loadSpeciesAbundancesFiltered.R")

# List all possible sample pairs.
# Do not include blacklisted samples.
samplePairsAll <- samplesRaw %>%
  dplyr::select(sample) %>% dplyr::rename(sample1=sample) %>%
  filter(!sample1 %in% sampleBlacklist) %>%
  mutate(sample2=sample1) %>% expand(sample1, sample2)

# Filter the list of possible sample pairs to include only pairs of samples
# from different subjects in subject X at the initial timepoint.
samplePairsInitialX <- samplePairsAll %>%
  filter(sample1 %in% samplesInitial, sample2 %in% samplesInitial,
         substr(sample1,1,3) %in% subjectsX, substr(sample2,1,3) %in% subjectsX)

# Write a function that, for a given pair of samples, calculates the number of species
# shared above the limit of detection and the percent relative abundance that the species make up.
options(dplyr.summarise.inform = FALSE)
compareSpecies <- function(sampleAbundance, sampleAnnotation){
  # Extract the species abundances for both samples.
  speciesAbundancesAbundance <- speciesAbundances %>%
    filter(sample==sampleAbundance)
  speciesAbundancesAnnotation <- speciesAbundances %>%
    filter(sample==sampleAnnotation)
  # Extract that list of species present above the limit of detection for the annotation sample.
  speciesAnnotation <- (speciesAbundancesAnnotation %>%
                          filter(relative_abundance>LIMITOFDETECTION))$species_id
  # Annotate the species present above the limit of detection in the abundance sample
  # with their presence or absence in the annotation sample.
  speciesAbundancesAbundance <- speciesAbundancesAbundance %>%
    filter(relative_abundance>LIMITOFDETECTION) %>%
    mutate(speciesShared=(species_id %in% speciesAnnotation), annotationSample=sampleAnnotation) %>%
    dplyr::select(sample, annotationSample, species_id, relative_abundance, speciesShared)
  # Summarize the number of species shared and their relative abundance.
  speciesAbundancesAbundanceSummary <- speciesAbundancesAbundance %>%
    mutate(totalSpecies=n()) %>%
    group_by(sample, annotationSample, totalSpecies, speciesShared) %>%
    summarize(numSpeciesShared=n(), relAbundanceSpeciesShared=sum(relative_abundance))
  return(speciesAbundancesAbundanceSummary)
}

# Compare the number and proportion of species shared for all initial samples
# from study arm X.
iSpeciesComparisons=0
speciesComparisonsInitial <- foreach(sampleAbundance=samplePairsInitialX$sample1,
                                     sampleAnnotation=samplePairsInitialX$sample2, .combine="rbind") %do%
  {
    iSpeciesComparisons <- iSpeciesComparisons+1
    if(iSpeciesComparisons %% 100 == 0){
      print(paste0(iSpeciesComparisons, " species comparisons performed"))
    }
    return(compareSpecies(sampleAbundance, sampleAnnotation))
  }
# Annotate species comparisons as within-subject, within-household, or between-household
# and plot the distribution of values.
speciesComparisonsInitial <- speciesComparisonsInitial %>%
  mutate(hh1=substr(sample,1,2), hh2=substr(annotationSample,1,2)) %>%
  mutate(type=ifelse(sample==annotationSample, "sameSubject",
                     ifelse(hh1==hh2, "sameHousehold","diffHousehold")))

# Export the species comparisons.
write.table(speciesComparisonsInitial, paste0(OUTDIR, "speciesComparisons-initial-X.txt"),
            quote=FALSE, row.names=FALSE)