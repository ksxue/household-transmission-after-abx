# This script takes pairs of samples at key timepoints from cohabiting subjects,
# and annotates each species in the abundance sample based on whether it is
# present in the annotation sample and whether it is shared with the abundance sample.

library(tidyverse)
library(foreach)

# Import metadata on the samples and subjects.
source("workflow/analysis/background/background.R")
# Set the output directory.
OUTDIR <- "workflow/analysis/integrateSpeciesAbundancesSharedStrains/out/"
# Set the limit of detection for species detected through metagenomics.
LIMITOFDETECTION <- 0.001

# Import the list of species abundances.
source("workflow/analysis/generateSpeciesAbundances/loadSpeciesAbundancesFiltered.R")
# Import the list of strain sharing calls aggregated from fixed differences and strain fishing
# for the initial samples from study arm X.
dataStrainSharing <- read.table("workflow/analysis/compareSharedStrainCalls/out/strainSharingCalls-consolidated-sameSubject.txt",
                                header=TRUE, stringsAsFactors = FALSE)

# Import the list of strain fishing calculations.
# These metrics help quantify the proportion of a population that's shared between two samples.
dataStrainFishing <- read.table("workflow/analysis/identifySharedStrains-strainFishing/out/strainFishing-sameSubject.txt.gz",
                                header=TRUE, stringsAsFactors = FALSE)
# Import the frequencies of bait strains.
# These metrics are also helpful for quantifying the proportion of a population that's shared between two samples.
dataBaitStrains <- read.table("workflow/analysis/identifySharedStrains-strainFishing/out/baitStrains.txt",
                              header=TRUE, stringsAsFactors = FALSE)

# Write a function to annotate species abundances with strain sharing information.
annotateSpeciesAbundancesSharedStrains <- function(sampleAbundance, sampleAnnotation, limitOfDetection){

  # Annotate each species in the abundance sample based on whether it is shared with
  # the annotation sample at both the species and strain level.
  # Extract the alphabetically earlier sample for use with the strain sharing dataframe.
  sample1value <- ifelse(sampleAbundance<sampleAnnotation, sampleAbundance, sampleAnnotation)
  # Extract the species abundances for both samples.
  speciesAbundancesAbundance <- speciesAbundances %>%
    filter(sample==sampleAbundance)
  speciesAbundancesAnnotation <- speciesAbundances %>%
    filter(sample==sampleAnnotation)
  # Extract that list of species present above the limit of detection for the annotation sample.
  speciesAnnotation <- (speciesAbundancesAnnotation %>%
                          filter(relative_abundance>limitOfDetection))$species_id
  # Extract the set of strain sharing calls for this pair of timepoints.
  strainSharingAnnotations <- dataStrainSharing %>%
    filter((sample1==sampleAbundance & sample2==sampleAnnotation) | 
             (sample1==sampleAnnotation & sample2==sampleAbundance)) %>%
    dplyr::rename(species_id=species) %>%
    dplyr::select(species_id, sample1, sample2,
                  consensusConservative, consensusPermissive, 
                  consensusFixedDiffs, consensusStrainFishing, consensusNA)
  # Change the orientation of the strain sharing annotations to match the "orientation"
  # of the abundance-annotation pair.
  if(sampleAnnotation==sample1value){
    strainSharingAnnotations <- strainSharingAnnotations %>%
      dplyr::mutate(sample=sample2, sample2=sample1, sample1=sample) %>%
      dplyr::select(-sample)
  }
  # Merge the species abundances in the abundance sample with the strain sharing annotations.
  speciesAbundancesAbundanceAnnotated <- 
    left_join(speciesAbundancesAbundance, 
              strainSharingAnnotations %>% 
                dplyr::rename(sample=sample1), 
              by=c("sample", "species_id"))
  # Parse the strain sharing annotations and also add annotations about species presence/absence.
  speciesAbundancesAbundanceAnnotated <- speciesAbundancesAbundanceAnnotated %>%
    mutate(annotationSample=sampleAnnotation,
           speciesPresentAbundanceSample=(relative_abundance>0),
           speciesPresentAnnotationSample=(species_id %in% speciesAnnotation),
           speciesSharing=
             ifelse(speciesPresentAbundanceSample & speciesPresentAnnotationSample, "speciesShared",
                    ifelse(speciesPresentAbundanceSample & (!speciesPresentAnnotationSample), "speciesNotShared", "speciesAbsent")),
           strainSharing=ifelse(is.na(consensusNA), "strainsUnknown",
                                ifelse(consensusNA, "strainsShared", "strainsNotShared"))) %>%
    dplyr::select(sample, annotationSample, species_id, relative_abundance, 
                  speciesPresentAbundanceSample, speciesPresentAnnotationSample, speciesSharing,
                  consensusNA, strainSharing) %>%
    mutate(annotation=ifelse(strainSharing=="strainsUnknown",
                             ifelse(speciesSharing=="speciesShared", "strainsUnknown", "speciesNotShared"), strainSharing))
  
  # For species that share strains, calculate the proportion of the overall population of the species
  # that the shared strain represents.
  # In cases of strain sharing, extract the proportion of each population in the sampleAbundance
  # that is shared with the sampleAnnotation based on the strain fishing calculations.
  strainFishingPair <- dataStrainFishing %>%
    filter((bait==sampleAnnotation & sample==sampleAbundance) |
             (bait==sampleAbundance & sample==sampleAnnotation)) %>%
    filter(shared)
  # Append the frequencies of the bait samples to the strain fishing data about these two samples.
  strainFishingPair <- left_join(strainFishingPair,
                                 dataBaitStrains %>% 
                                   dplyr::rename(baitMeanSNPfreq=meanSNPfreq) %>%
                                   dplyr::select(species, bait, baitMeanSNPfreq),
                                 by=c("species","bait"))
  # Extract the frequencies of each population that are made up of the shared strains.
  strainFishingPairFrequencies <- strainFishingPair %>%
    mutate(sharedStrainFreq=ifelse(bait==sampleAbundance, baitMeanSNPfreq, medianSNPfreq)) %>%
    dplyr::select(species, bait, sample, baitMeanSNPfreq, medianSNPfreq, sharedStrainFreq)
  # If both directions of the strain fishing comparison have been calculated,
  # then use the shared strain frequency calculated using the sampleAnnotation as the bait.
  strainFishingPairFrequencies <- strainFishingPairFrequencies %>%
    group_by(species) %>% 
    mutate(annotationSample=sampleAnnotation,
           numComparisons=n(),
           mainAnnotation=ifelse(numComparisons==1,TRUE,(bait==annotationSample))) %>%
    filter(mainAnnotation) %>% dplyr::select(species, sharedStrainFreq)
  
  # Append the frequencies of the shared strains to the annotated relative abundance dataframe.
  speciesAbundancesAbundanceAnnotated <- 
    left_join(speciesAbundancesAbundanceAnnotated, 
              strainFishingPairFrequencies %>% dplyr::rename(species_id=species),
              by=c("species_id"))
  
  return(speciesAbundancesAbundanceAnnotated)
}

# List all possible sample pairs that involve the initial and final timepoint
# from the same subject.
# Do not include blacklisted samples.
samplePairsAll <- samplesRaw %>%
  dplyr::select(sample) %>% dplyr::rename(sample1=sample) %>%
  filter(!sample1 %in% sampleBlacklist) %>%
  mutate(sample2=sample1) %>% expand(sample1, sample2)
# Retain only sample pairs from the initial timepoint,
# including both cohabiting and noncohabiting pairs.
samplePairsInitialFinal<- samplePairsAll %>%
  filter(((sample1 %in% c(samplesInitial) & sample2 %in% c(samplesFinal)) |
            (sample2 %in% c(samplesInitial) & sample1 %in% c(samplesFinal))),
         substr(sample1,1,1)=="X", substr(sample2,1,1)=="X",
         substr(sample1,1,3)==substr(sample2,1,3))


# Annotate the species abundances for each within-household pair of samples.
iSamplePairs=0
speciesAbundancesAnnotatedInitialFinal <- 
  foreach(isample1=samplePairsInitialFinal$sample1,
          isample2=samplePairsInitialFinal$sample2, .combine="rbind") %do% {
            # Increment counter.
            iSamplePairs <- iSamplePairs+1
            if(iSamplePairs %% 100 == 0){
              print(paste0(iSamplePairs, " sample pairs analyzed"))
            }
            # Return output of annotation function.
            return(annotateSpeciesAbundancesSharedStrains(isample1, isample2, LIMITOFDETECTION))
          }

# Export the annotated species abundances for within-household samples.
# Remove rows corresponding to species with relative abundance 0 to save space.
write.table(speciesAbundancesAnnotatedInitialFinal %>% filter(relative_abundance>0),
            gzfile(paste0(OUTDIR, "annotatedSpeciesAbundances-sameSubject-initialFinal.txt.gz")),
            row.names=FALSE, quote=FALSE)

# For each sample pair, summarize the number of species and the percent relative abundance
# in each of the annotation categories.
speciesAbundancesAnnotatedInitialFinalSummary <-
  speciesAbundancesAnnotatedInitialFinal %>%
  filter(relative_abundance>0) %>%
  group_by(sample, annotationSample, annotation) %>%
  summarize(numSpecies=n(), totRelAbundance=sum(relative_abundance))
# Export the summary of the relative abundances in each annotation class.
write.table(speciesAbundancesAnnotatedInitialFinalSummary,
            paste0(OUTDIR, "annotatedSpeciesAbundances-sameSubject-initialFinal-summary.txt"),
            row.names=FALSE, quote=FALSE)


