library(tidyverse)
library(data.table)

# Import plot helpers.
source("workflow/analysis/plotDefaults.R")
OUTDIR <- "workflow/analysis/calculateHMPspeciesPrevalenceAbundance/out/"

# Import the sample metadata.
source("workflow/analysis/background/background.R")

# Import the HMP species abundances from Good and Rosenberg,
# which were calculated using MIDAS species profiling from the HMP.
# This dataframe includes 932 samples.
speciesAbundancesHMP <- fread("data/HMP/coverage.txt",
                              header=TRUE, stringsAsFactors = FALSE, sep="\t", data.table = FALSE)
# Import the metadata for the HMP samples.
HMPmetadata <- fread("data/HMP/HMP1-2_ids_order.txt",
                     stringsAsFactors = FALSE, data.table=FALSE)
# There are 553 samples from 249 subjects.
n_distinct(HMPmetadata$subject_id)
n_distinct(HMPmetadata$sample_id)

# Tidy the full HMP data.
speciesAbundancesHMPtidy <- speciesAbundancesHMP %>%
  pivot_longer(-species_id, names_to="sample_id", values_to="coverage")
# Append the metadata and filter to include only samples from the HMP.
speciesAbundancesHMPtidy <- 
  left_join(speciesAbundancesHMPtidy %>% filter(sample_id %in% HMPmetadata$sample_id) %>%
              mutate(sample_id=as.integer(sample_id)), 
            HMPmetadata, by="sample_id")

# For each sample, calculative the relative abundance of each species.
# There are 391 samples that are recognized as part of the HMP.
speciesAbundancesHMPtidy <- speciesAbundancesHMPtidy %>%
  group_by(sample_id) %>%
  mutate(relative_abundance=coverage/sum(coverage))
# Retain only a single sample from each subject.
# Keep the one with the smaller sample_id.
# This results in 229 samples from distinct subjects in the HMP.
speciesAbundancesHMPPerSubject <- speciesAbundancesHMPtidy %>%
  ungroup() %>% group_by(subject_id) %>%
  filter(sample_id==min(sample_id))

# Calculate the prevalence and abundance of each species in the HMP cohort.
# Define the prevalence as the number of subjects in which the species is detected
# at any frequency, and the abundance as the median relative abundance in
# subjects that have the species at a non-zero frequency.
# Remove species that have 0 prevalence in this dataset.
speciesPrevalenceAbundanceHMP <- speciesAbundancesHMPPerSubject %>%
  ungroup() %>% group_by(species_id) %>%
  summarize(prevalence=sum(relative_abundance>0)/n(), 
            medianAbundance=median(relative_abundance[relative_abundance>0]),
            meanAbundance=mean(relative_abundance[relative_abundance>0]),
            maxAbundance=max(relative_abundance)) %>%
  filter(prevalence>0)
# Export the prevalence and abundance.
write.table(speciesPrevalenceAbundanceHMP, paste0(OUTDIR, "HMP-speciesPrevalenceAbundance.txt"),
            quote=FALSE, row.names=FALSE)

# Export the species with high prevalence and abundance.
write.table(speciesPrevalenceAbundanceHMP %>% filter(prevalence>0.05, medianAbundance>1e-4) %>%
              arrange(desc(prevalence)),
            paste0(OUTDIR, "HMP-coreSpecies.txt"))
