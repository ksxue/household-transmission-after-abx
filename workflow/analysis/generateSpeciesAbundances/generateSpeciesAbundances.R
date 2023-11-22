library(data.table)
library(tidyverse)

# Import background file with all sample and study information.
source("workflow/analysis/background/background.R")

# List input file with species profiles for each sample.
INPUT.SPECIESPROFILE <- "workflow/out/midasOutput/species/species_profile_all.txt"
# List input file with species taxonomy.
INPUT.SPECIESTAXONOMY <- "workflow/out/midasOutput/database/species_taxonomy.txt"
# Output directory
OUTDIR <- "workflow/analysis/generateSpeciesAbundances/out/"

# Import species abundance and taxonomy data.
# Import species relative abundances.
data <- fread(INPUT.SPECIESPROFILE, 
              header=TRUE, stringsAsFactors = FALSE, data.table=FALSE)
# Complete the table of species abundances.
# Populate missing observations with zero.
data <- data %>%
  complete(sample, species_id, 
           fill=list(relative_abundance=0, count_reads=0, coverage=0))
# Parse the sample names.
data <- data %>%
  mutate(sample=substr(sample,29,35),
         subject=substr(sample,1,3),
         timepoint=as.numeric(substr(sample,5,7)))
# Exclude samples from undefined subjects.
data <- data %>%
  filter(sample %in% samplesAll)
# Calculate the total number of reads mapped to single-copy genes
# for each sample.
data <- data %>%
  group_by(sample) %>%
  mutate(totalReads=sum(count_reads))

# Import species taxonomies from MIDAS genome database.
speciesTaxons <- fread(INPUT.SPECIESTAXONOMY,
                       header=TRUE, stringsAsFactors = FALSE, data.table=FALSE)
speciesTaxons <- speciesTaxons %>% 
  dplyr::select(-genome_id,-genome_name,-taxon_id,
                -taxon_lineage_ids,-taxon_lineage_names)
# Shorten some long family names.
speciesTaxons <- speciesTaxons %>%
  mutate(family=gsub("\\..*","",family))
# Add species taxonomy to relative abundances.
data <- left_join(data, 
                  speciesTaxons,
                  by=c("species_id"))


# Summarize the total coverage for each sample
# in terms of reads mapped to single-copy genes.
dataCoverage <- data %>%
  group_by(sample) %>%
  summarize(totalReads=sum(count_reads))
# Export species sequencing coverage.
write.table(dataCoverage, paste0(OUTDIR,"speciesCoverage.txt"),
            quote=FALSE, row.names = FALSE)

# Export species abundances. Remove species with zero counts to save space.
write.table(data %>% filter(count_reads>0), paste0(OUTDIR,"speciesAbundances.txt"),
            quote=TRUE, row.names = FALSE)
