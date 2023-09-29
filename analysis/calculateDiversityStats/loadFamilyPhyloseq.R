# Imports the phyloseq object generated from the species abundance table.

library(tidyverse)
library(phyloseq)

# Import the phyloseq object.
ps_family_filtSample <- readRDS("workflow/analysis/calculateDiversityStats/out/families-filtSample.ps")
ps_family_filtSample_filtAbundance <- 
  readRDS("workflow/analysis/calculateDiversityStats/out/families-filtSample-filtAbundance.ps")
