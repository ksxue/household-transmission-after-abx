# This code loads species abundances and prepares them into a standard format
# for plotting and analyses.
# It also filters out samples that are blacklisted due to contamination
# or suspected sample swaps.

library(data.table)
library(tidyverse)

# This code loads species abundances and prepares them into a standard format
# for plotting and analyses.
# Set the limit of detection for filtering species abundances.
LIMITOFDETECTION <- 1e-3

# Import raw species abundances.
speciesAbundances <- fread("workflow/analysis/generateSpeciesAbundances/out/speciesAbundances.txt",
                           header = TRUE, data.table = FALSE)
# Fill out the species abundances such that the dataframe is no longer ragged.
speciesAbundances <- speciesAbundances %>%
  complete(sample, species_id, 
           fill=list(count_reads=0, coverage=0, relative_abundance=0, totalReads=0))

# Import species taxonomies from MIDAS genome database.
# Restore these taxonomies for the species that were added using "complete."
speciesTaxons <- fread("workflow/out/midasOutput/database/species_taxonomy.txt",
                       header=TRUE, stringsAsFactors = FALSE, data.table=FALSE)
speciesTaxons <- speciesTaxons %>% 
  dplyr::select(-genome_id,-genome_name,-taxon_id,
                -taxon_lineage_ids,-taxon_lineage_names)
# Shorten some long family names.
speciesTaxons <- speciesTaxons %>%
  mutate(family=gsub("\\..*","",family))
# Add species taxonomy to relative abundances.
speciesAbundances <- left_join(speciesAbundances %>% 
                                 dplyr::select(sample, species_id, count_reads,
                                               coverage, relative_abundance, totalReads), 
                               speciesTaxons,
                               by=c("species_id"))

# Parse sample name information.
speciesAbundances <- speciesAbundances %>%
  mutate(subject=substr(sample,1,3),
         timepoint=as.numeric(substr(sample,5,7)),
         hh=substr(sample,1,2))

# Generate a species abundance dataframe that lists only species above 0.1%.
# This will help control for differences in sample sequencing depth.
speciesAbundancesFiltAbundance <- speciesAbundances %>%
  filter(relative_abundance>LIMITOFDETECTION) %>%
  complete(sample, species_id, 
           fill=list(relative_abundance=0, count_reads=0, coverage=0))
# Add species taxonomy to relative abundances.
speciesAbundancesFiltAbundance <- left_join(speciesAbundancesFiltAbundance %>% 
                                 dplyr::select(sample, species_id, count_reads,
                                               coverage, relative_abundance, totalReads), 
                               speciesTaxons,
                               by=c("species_id"))
# Parse sample name information.
speciesAbundancesFiltAbundance <- speciesAbundancesFiltAbundance %>%
  mutate(subject=substr(sample,1,3),
         timepoint=as.numeric(substr(sample,5,7)),
         hh=substr(sample,1,2))

# Load list of blacklisted samples.
sampleBlacklist <- (read.table("workflow/analysis/background/out/sampleBlacklist.txt",
                               header=TRUE, stringsAsFactors = FALSE))$sample

# Filter out samples that are on the sample blacklist.
speciesAbundances <- speciesAbundances %>%
  filter(!(sample %in% sampleBlacklist))
speciesAbundancesFiltAbundance <- speciesAbundancesFiltAbundance %>%
  filter(!(sample %in% sampleBlacklist))

# Identify samples that have lower than 10^3.5 reads.
lowCovSamples <- speciesAbundances %>% 
  dplyr::select(sample, totalReads) %>% unique() %>% 
  arrange(totalReads) %>% filter(totalReads>0, totalReads<10^3.5)
# Filter out low-coverage samples.
speciesAbundances <- speciesAbundances %>%
  group_by(sample) %>% mutate(totalReads=sum(count_reads)) %>%
  ungroup() %>% filter(totalReads>10^3.5)
speciesAbundancesFiltAbundance <- speciesAbundancesFiltAbundance %>%
  group_by(sample) %>% mutate(totalReads=sum(count_reads)) %>%
  ungroup() %>% filter(totalReads>10^3.5)

