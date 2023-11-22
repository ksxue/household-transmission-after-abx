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
generatePhyloseq <- function(xSpeciesAbundances){
  # Generate OTU table from species abundances.
  dataOTU <- xSpeciesAbundances %>% ungroup() %>%
    dplyr::select(sample, species_id, count_reads) %>%
    spread(sample, count_reads) %>% arrange(species_id)
  rownames(dataOTU) <- dataOTU$species_id
  speciesList <- dataOTU$species_id
  dataOTU <- dataOTU %>% dplyr::select(-species_id)
  dataOTU <- as.matrix(dataOTU)
  rownames(dataOTU) <- speciesList
  
  # Generate taxonomy table from MIDAS genome database taxonomy.
  dataTax <- xSpeciesAbundances %>% ungroup() %>%
    dplyr::select(species_id, kingdom, phylum, class, order, family, genus, species) %>%
    dplyr::rename(Domain=kingdom, Phylum=phylum, Class=class,
                  Order=order, Family=family, Genus=genus, Species=species) %>%
    unique() %>% arrange(species_id)
  rownames(dataTax) <- dataTax$species_id
  dataTax <- dataTax %>% dplyr::select(-species_id)
  dataTax <- as.matrix(dataTax)
  rownames(dataTax) <- speciesList
  
  # Combine objects into a phyloseq object.
  ps <- phyloseq(otu_table(dataOTU, taxa_are_rows = TRUE),
                 tax_table(dataTax))
  
  return(ps)
}


# Export phyloseq object.
saveRDS(generatePhyloseq(speciesAbundances), 
        paste0(OUTDIR, "species-filtSample.ps"))
saveRDS(generatePhyloseq(speciesAbundancesFiltAbundance), 
        paste0(OUTDIR, "species-filtSample-filtAbundance.ps"))
