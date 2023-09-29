# This script takes the filtered dataframe of fixed differences
# between pairs of samples of each high-coverage species.
# It uses the between-household distribution of fixed differences
# to set empirical thresholds for what constitutes a shared strain,
# and it identifies strains within households 
# that are more similar than this threshold.

library(tidyverse)
library(data.table)

# Import background information on study design.
source("workflow/analysis/background/background.R")
# Import the filtered fixed differences dataframe.
source("workflow/analysis/identifySharedStrains-fixedDiffs/loadFixedDiffsFiltered.R")
# Set an output directory.
OUTDIR <- "workflow/analysis/identifySharedStrains-fixedDiffs/out/"

# Filter out samples that have a median sequencing coverage below 5,
# since they are likely to have many fixed differences due to noise.
dataFixedDiffs <- dataFixedDiffs %>%
  filter(sample1cov>=5, sample2cov>=5)

# For each species, calculate the minimum number of fixed differences 
# between pairs of strains from subjects in different households.
dataFixedDiffs <- dataFixedDiffs %>%
  group_by(species) %>%
  mutate(btwnThreshold=ifelse(sum(type=="diffHousehold")==0, NA, 
                              round(quantile(fixedDiffs[type=="diffHousehold"], probs=c(0.01)))))
# Summarize the between-household thresholds for number of fixed differences between pairs of samples.
btwnThresholds <- dataFixedDiffs %>%
  dplyr::select(species, btwnThreshold) %>%
  unique()

# Annotate the list of fixed differences.
# Identify pairs of populations that harbor shared strains.
dataFixedDiffsShared <- dataFixedDiffs %>%
  filter(type=="sameHousehold") %>%
  mutate(shared=(fixedDiffs<btwnThreshold))

# Export strain sharing calls within subjects, retaining all metadata for downstream analyses.
write.table(dataFixedDiffs %>% filter(type=="sameSubject") %>%
              mutate(shared=(fixedDiffs<btwnThreshold)),
            gzfile(paste0(OUTDIR, "fixedDiffs-sameSubject.txt.gz")),
            row.names=FALSE, quote=FALSE)

# Export strain sharing calls, retaining all metadata for downstream analyses.
write.table(dataFixedDiffsShared,
            gzfile(paste0(OUTDIR, "fixedDiffs-sharingAnnotated.txt.gz")),
            row.names=FALSE, quote=FALSE)

# # Export the fixed differences data for all pairs of samples with coverage >5.
# write.table(dataFixedDiffs,
#             gzfile(paste0(OUTDIR, "fixedDiffs-allSamples.txt.gz")),
#             row.names=FALSE, quote=FALSE)

# Export the strain sharing calls for all samples from the initial timepoint.
write.table(dataFixedDiffs %>%
              filter(sample1 %in% samplesInitial, sample2 %in% samplesInitial,
                     subject1 %in% subjectsX, subject2 %in% subjectsX, sample1!=sample2) %>%
              mutate(shared=(fixedDiffs<btwnThreshold)),
            gzfile(paste0(OUTDIR, "fixedDiffs-initialX.txt.gz")),
            row.names=FALSE, quote=FALSE)

# Import the genome length for each species to calculate the % fixed differences.
# Import species taxonomies from MIDAS genome database.
# Restore these taxonomies for the species that were added using "complete."
speciesTaxons <- fread("workflow/out/midasOutput/database/species_taxonomy.txt",
                       header=TRUE, stringsAsFactors = FALSE, data.table=FALSE)
genomeInfo <- fread("workflow/out/midasOutput/database/genome_info.txt",
                    header=TRUE, stringsAsFactors = FALSE, data.table=FALSE)
btwnThresholdsAnnotated <- left_join(
  left_join(btwnThresholds %>% dplyr::rename(species_id=species),
            speciesTaxons %>% dplyr::select(species_id, genome_id), by="species_id"),
  genomeInfo %>% dplyr::select(genome_id, length), by="genome_id")
btwnThresholdsAnnotated <- btwnThresholdsAnnotated %>%
  mutate(pctGenomeLength=btwnThreshold/length)

# Export the between-household thresholds for calling strains as shared.
write.table(btwnThresholdsAnnotated,
            (paste0(OUTDIR, "fixedDiffs-btwnHouseholdThresholds.txt")),
            row.names=FALSE, quote=FALSE)


