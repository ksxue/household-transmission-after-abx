# Calculate multiple alpha-diversity statistics
# for each of the sequenced samples.
# Omit samples that have previously been black-listed.
# Export a table listing multiple alpha-diversity stats for each sample.

library(tidyverse)
library(phyloseq)

# Import the phyloseq object.
source("workflow/analysis/calculateDiversityStats/loadSpeciesPhyloseq.R")
# Import the sample metadata.
source("workflow/analysis/background/background.R")
# Set the output directory.
OUTDIR <- "workflow/analysis/calculateDiversityStats/out/"

# Calculate alpha diversity.
calculateAlphaDiversity <- function(ps) {

  # Use the estimate_richness function to calculate
  # an array of alpha-diversity statistics.
  alphaRaw <- estimate_richness(ps, split=TRUE)
  # Add the sample names.
  alphaRaw$sample <- rownames(alphaRaw)
  # Calculate the Shannon effective number of species.
  alphaRaw <- alphaRaw %>%
    mutate(ShannonEffectiveSpecies=exp(Shannon))
  # Tidy the dataframe.
  alpha <- alphaRaw %>%
    mutate(sample=gsub("\\.","-",sample)) %>%
    pivot_longer(-sample, names_to="alphaStat", values_to="value")
  
}

# Calculate alpha diversity for both the full species dataset
# and the abundance-filtered species dataset.
speciesAlpha <- 
  rbind(calculateAlphaDiversity(ps_filtSample) %>%
          mutate(type="filtSample"),
        calculateAlphaDiversity(ps_filtSample_filtAbundance) %>%
          mutate(type="filtSample_filtAbundance"))

# Export the alpha diversity data.
write.table(speciesAlpha, paste0(OUTDIR, "speciesAlpha.txt"),
            quote=FALSE, row.names=FALSE)
