# This code loads filtered species abundances from MIDAS
# and converts them into phyloseq format for downstream analyses.
# Samples that have been blacklisted are removed from this analysis.

library(data.table)
library(tidyverse)
library(phyloseq)

# Import dataframe of species abundances.
source("workflow/analysis/generateSpeciesAbundances/loadSpeciesAbundancesFiltered.R")
# Set path of output file.
OUTDIR <- "workflow/analysis/calculateDiversityStats/out/"
LIMITOFDETECTION <- 1e-3


# Write a function to generate a phyloseq object from variants of
# the species abundance dataframe.
# Summarize taxa at the level of bacterial families.
generatePhyloseq <- function(xSpeciesAbundances){
  # Generate OTU table from species abundances.
  dataOTU <- xSpeciesAbundances %>% ungroup() %>%
    group_by(sample, family) %>%
    summarize(count_reads=sum(count_reads)) %>%
    dplyr::select(sample, family, count_reads) %>%
    spread(sample, count_reads) %>% arrange(family)
  rownames(dataOTU) <- dataOTU$family
  speciesList <- dataOTU$family
  dataOTU <- dataOTU %>% dplyr::select(-family)
  dataOTU <- as.matrix(dataOTU)
  rownames(dataOTU) <- speciesList
  
  # # Generate taxonomy table from MIDAS genome database taxonomy.
  # dataTax <- xSpeciesAbundances %>% ungroup() %>%
  #   dplyr::select(kingdom, phylum, class, order, family) %>%
  #   dplyr::rename(Domain=kingdom, Phylum=phylum, Class=class,
  #                 Order=order, Family=family) %>%
  #   unique() %>% arrange(Family)
  # rownames(dataTax) <- dataTax$Family
  # dataTax <- dataTax %>% dplyr::select(-family)
  # dataTax <- as.matrix(dataTax)
  # rownames(dataTax) <- speciesList

  # Combine objects into a phyloseq object.
  ps <- phyloseq(otu_table(dataOTU, taxa_are_rows = TRUE))
  
  return(ps)
}


# Export phyloseq object.
saveRDS(generatePhyloseq(speciesAbundances), 
        paste0(OUTDIR, "families-filtSample.ps"))
saveRDS(generatePhyloseq(speciesAbundancesFiltAbundance), 
        paste0(OUTDIR, "families-filtSample-filtAbundance.ps"))
