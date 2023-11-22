library(tidyverse)
library(data.table)

# If TEST is true, then run on default files
# rather than snakemake inputs.
TEST <- FALSE

# Import species info from MIDAS genome database.
speciesInfo <- read.table(ifelse(TEST,"/home/kxue/bin/MIDASdb/midas_db_v1.2/species_info.txt",
                                 snakemake@input[["speciesInfo"]]),
                          header=TRUE, stringsAsFactors = FALSE, fill=TRUE)

# Import genome taxonomies from MIDAS genome database.
genomeTaxons <- fread(ifelse(TEST,"/home/kxue/bin/MIDASdb/midas_db_v1.2/genome_taxonomy.txt",
                              snakemake@input[["genomeTaxonomy"]]),
                       header=TRUE, data.table=FALSE)

# Retain only genomes that are maintained as the representative
# genome for their species.
speciesTaxons <- genomeTaxons %>%
  filter(genome_id %in% speciesInfo$rep_genome)

# Merge the species names and taxa tables.
speciesTaxons <- left_join(speciesInfo %>% dplyr::select(species_id, rep_genome) %>%
                             dplyr::rename(genome_id=rep_genome),
                           speciesTaxons, by=c("genome_id"))

# Export the taxonomies of the species in the MIDAS database.
write.table(speciesTaxons,
            ifelse(TEST,"workflow/out/midasOutput/database/species_taxonomy.txt",
                   snakemake@output[["speciesTaxonomy"]]),
            quote=FALSE, row.names=FALSE, sep="\t")
