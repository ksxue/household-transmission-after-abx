# This script reads in the list of fixed differences
# between pairs of samples, and it identifies samples
# that have evidence of cross-contamination
# (that is, pairs of samples from different households
# that have very low numbers of fixed differences).

library(tidyverse)
library(data.table)

# Import the list of fixed differences between pairs of samples.
source("workflow/analysis/identifyContaminatedSamples/loadFixedDiffs.R")
# Set output directory.
OUTDIR <- "workflow/analysis/identifyContaminatedSamples/out/"

# Filter out samples that have a median sequencing coverage below 5,
# since they are likely to have many fixed differences due to noise.
dataFixedDiffsFiltered <- dataFixedDiffs %>%
  filter(sample1cov>=5, sample2cov>=5)

# Identify pairs of between-household samples 
# that have identical or near-identical fixed differences.
# These samples may be indications of contamination.
suspiciousPairsOfSamples <- dataFixedDiffsFiltered %>%
  filter(species %in% highCovSpecies, type=="diffHousehold",
         fixedDiffs<3)
# Export the list of suspicious samples.
write.table(suspiciousPairsOfSamples, paste0(OUTDIR, "suspiciousPairs.txt"),
            quote=FALSE, row.names=FALSE)

# Filter out suspicious sample pairs for which only one species is suspicious.
# We expect that true contaminants would be flagged as suspicious pairs
# across multiple species.
suspiciousPairsOfSamplesFiltered <- suspiciousPairsOfSamples %>%
  group_by(sample1, sample2) %>% mutate(numSpecies=n()) %>%
  filter(numSpecies>1)

# Using the suspicious sample pairs, identify suspicious samples.
if(nrow(suspiciousPairsOfSamplesFiltered)>0){
  # Take the list of suspicious sample pairs.
  # Identify the most suspicious sample and remove it from the list.
  # Note its impact in terms of how many suspicious sample pairs are removed
  # by removing this sample.
  # Regenerate the list of suspicious sample pairs
  # Repeat this process until no suspicious sample pairs remain.
  tempSuspiciousPairsOfSamples <- suspiciousPairsOfSamplesFiltered
  suspiciousSamplesIterative <- c()
  suspiciousSamplesIterativeImpact <- c()
  while(nrow(tempSuspiciousPairsOfSamples)>0){
    # Identify the most suspicious sample of these pairs.
    # Note that top_n will return multiple samples if there is a tie.
    tempMostSuspiciousSamples <- (tibble(sample=c(tempSuspiciousPairsOfSamples$sample1, 
                                                  tempSuspiciousPairsOfSamples$sample2)) %>%
                                    group_by(sample) %>% summarize(numSuspiciousPairs=n()) %>% 
                                    top_n(1, numSuspiciousPairs))$sample
    # To break the tie, designate the sample that comes alphabetically last
    # as the most suspicious sample. This reflects the fact that samples from
    # study arms Y and Z are more likely to have low DNA yield.
    tempMostSuspiciousSample <- (sort(tempMostSuspiciousSamples))[length(tempMostSuspiciousSamples)]
    # Remove all suspicious pairs containing the suspicious sample from the list.
    newSuspiciousPairsOfSamples <- tempSuspiciousPairsOfSamples %>%
      filter(sample1 != tempMostSuspiciousSample,
             sample2 != tempMostSuspiciousSample)
    # Record the most suspicious sample and the number of suspicious pairs it accounts for.
    suspiciousSamplesIterative <- c(suspiciousSamplesIterative, 
                                    tempMostSuspiciousSample)
    suspiciousSamplesIterativeImpact <- 
      c(suspiciousSamplesIterativeImpact,
        nrow(tempSuspiciousPairsOfSamples)-nrow(newSuspiciousPairsOfSamples))
    # Begin the next iteration.
    tempSuspiciousPairsOfSamples <- newSuspiciousPairsOfSamples
  }
  # Store the list of samples inferred to be contaminated
  # and the number of suspect pairs accounted for by each sample.
  contaminatedSamples <- as.data.frame(cbind(suspiciousSamplesIterative,
                                             suspiciousSamplesIterativeImpact))
  colnames(contaminatedSamples) <- c("sample","numSuspiciousPairs")
  # Export the list of contaminated samples.
  write.table(contaminatedSamples, paste0(OUTDIR, "contaminatedSamples.txt"),
              quote=FALSE, row.names=FALSE)
  
  # For each sample in the list of contaminated samples, determine which subject
  # accounts for the majority of the cross-contamination.
  # Extract samples from the list of suspicious sample pairs
  # that involve the sample in question.
  suspiciousSamplesPatterns <- suspiciousPairsOfSamplesFiltered %>%
    filter(sample1 %in% contaminatedSamples$sample | sample2 %in% contaminatedSamples$sample) %>%
    mutate(contaminatedSample=ifelse(sample1 %in% contaminatedSamples$sample, sample1, sample2),
           contaminatingSample=ifelse(sample1 %in% contaminatedSamples$sample, sample2, sample1),
           contaminatingSubject=substr(contaminatingSample,1,3)) %>%
    group_by(contaminatedSample, contaminatingSubject) %>% 
    summarize(numSamples=n())
  write.table(suspiciousSamplesPatterns, paste0(OUTDIR, "contaminatedSamples-patterns.txt"),
              quote=FALSE, row.names=FALSE)
}

