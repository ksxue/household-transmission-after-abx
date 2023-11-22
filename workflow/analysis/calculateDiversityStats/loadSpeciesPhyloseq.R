# Imports the phyloseq object generated from the species abundance table.

library(tidyverse)
library(phyloseq)

# Import the phyloseq object.
ps_filtSample <- readRDS("workflow/analysis/calculateDiversityStats/out/species-filtSample.ps")
ps_filtSample_filtAbundance <- 
  readRDS("workflow/analysis/calculateDiversityStats/out/species-filtSample-filtAbundance.ps")
