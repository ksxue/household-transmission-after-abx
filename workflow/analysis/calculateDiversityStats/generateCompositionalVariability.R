# This script takes the species abundance dataframe, subsets on key timepoints,
# and converts it into a Q matrix.
# It then calculates normalized FST using the FSTRUCT package.
# Some of this code is borrowed from Maike Morrison's analyses.
library(tidyverse)
library(FSTruct)
library(foreach)

# Import the background information on each sample.
source("workflow/analysis/background/background.R")
# Import the dataframe containing species abundances.
source("workflow/analysis/generateSpeciesAbundances/loadSpeciesAbundancesFiltered.R")
# Set output directory.
OUTDIR <- "workflow/analysis/calculateDiversityStats/out/"

# Write a function to convert species abundances to a Q matrix
convert_to_Q <- function(df){
  pivot_wider(data = df, 
              id_cols = timepoint, 
              names_from = species_id, 
              values_from = relative_abundance, 
              values_fill = 0)
}

# Calculate FST (compositional variability) for each subject.
# Use all timepoints from the main study that are available for each subject.
FSTMainTimepoints <- foreach(iSubject=subjectsX, .combine="rbind") %do% {
  # Extract the species-abundance data corresponding to one subject.
  iSpeciesAbundances <- speciesAbundances %>%
    filter(subject==iSubject, sample %in% samplesXmain)
  if(nrow(iSpeciesAbundances)>0){
    # Convert the species abundances into a Q-matrix.
    iSubjectQ <- convert_to_Q(iSpeciesAbundances)
    # Calculate FST and other statistics.
    return(as.data.frame(Q_stat(iSubjectQ, n_distinct(iSpeciesAbundances$species_id))) %>%
             mutate(subject=iSubject))
  }
}


# To ensure a fair comparison between subjects with different numbers
# of sequenced samples, compare only the four key timepoints that were sequenced
# deeply for all subjects.
# Extract data representing the key timepoints from each subject.
speciesAbundancesKeyTimepoints <- speciesAbundances %>%
  filter(sample %in% samplesKeyTimepoints)


# Calculate the normalized FST (compositional variability) for each subject.
# This matrix is calculated using only the key timepoints.
FSTKeyTimepoints <- foreach(iSubject=subjectsX, .combine="rbind") %do% {
  # Extract the species-abundance data corresponding to one subject.
  iSpeciesAbundances <- speciesAbundancesKeyTimepoints %>%
    filter(subject==iSubject)
  if(nrow(iSpeciesAbundances)>0){
    # Convert the species abundances into a Q-matrix.
    iSubjectQ <- convert_to_Q(iSpeciesAbundances)
    # Calculate FST and other statistics.
    return(as.data.frame(Q_stat(iSubjectQ, n_distinct(iSpeciesAbundances$species_id))) %>%
             mutate(subject=iSubject))
  }
}

# Export the compositional variability statistics calculated across all timepoints.
write.table(FSTMainTimepoints, paste0(OUTDIR,"speciesCV-mainStudy.txt"), 
            quote=FALSE, row.names=FALSE)
write.table(FSTKeyTimepoints, paste0(OUTDIR,"speciesCV-keyTimepoints.txt"), 
            quote=FALSE, row.names=FALSE)
