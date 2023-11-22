# Calculate beta diversity using a range of diversity statistics.
# Use the sample-filtered phyloseq object that excludes blacklisted samples.
# Calculate these statistics at both the species and family levels.
# Export a table containing the beta-diversity statistics.

library(tidyverse)
library(phyloseq)
library(foreach)

# Load the species- and family-abundance phyloseq objects.
source("workflow/analysis/calculateDiversityStats/loadSpeciesPhyloseq.R")
source("workflow/analysis/calculateDiversityStats/loadFamilyPhyloseq.R")
# Import the background files.
source("workflow/analysis/background/background.R")
# Set the output directory.
OUTDIR <- "workflow/analysis/calculateDiversityStats/out/"

# List the distance-calculation methods to be used.
distMethods <- c("jsd","bray","jaccard")

# Write a function to calculate a distance matrix using the specified method
# and convert the data into tidy format.
calculateBeta <- function(ps, distMethod) {
  # Calculate the distance matrix using the specified method.
  betaRaw <- distance(ps, method=distMethod)
  
  # Convert distance matrix to a dataframe.
  beta <- as.matrix(betaRaw)
  beta <- as.data.frame(beta)
  beta$sample1 <- rownames(beta)
  # Tidy the dataframe.
  beta <- beta %>%
    pivot_longer(-sample1, names_to="sample2", values_to="value")
  beta <- beta %>%
    filter(sample1 != sample2) %>%
    mutate(method=distMethod)
}

# Calculate the distance matrix for all of the specified methods on the species abundances.
# Combine the distance matrices for all methods.
betaSpecies <- foreach(distMethod=distMethods, .combine = "rbind") %do% {
  print(distMethod)
  calculateBeta(ps_filtSample, distMethod)
}
# Export the distance matrix generated for all of the sample pairs
# using all of the specified methods.
write_delim(betaSpecies, paste0(OUTDIR, "speciesBeta.txt.gz"))

# Calculate the distance matrix for all of the specified methods on the family abundances.
# Combine the distance matrices for all methods.
betaFamily <- foreach(distMethod=distMethods, .combine = "rbind") %do% {
  print(distMethod)
  calculateBeta(ps_family_filtSample, distMethod)
}
# Export the distance matrix generated for all of the sample pairs
# using all of the specified methods.
write_delim(betaFamily, paste0(OUTDIR, "familiesBeta.txt.gz"))

