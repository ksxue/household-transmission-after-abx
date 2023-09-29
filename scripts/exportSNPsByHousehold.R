# This script takes the MIDAS output files, filters out sites with abnormal coverage,
# and calculates the number of fixed differences between each pair of samples.

library(data.table)
library(tidyverse)
library(cowplot)
library(foreach)
library(forcats)
library(R.utils)

# Import helper functions for parsing the MIDAS output files.
source("workflow/scripts/helperFunctions.R")

# If TEST is true, then run on default files
# rather than snakemake inputs.
TEST <- FALSE
if(TEST){
  setwd("Y:/oak/stanford/groups/relman/users/kxue/household-transmission-mgx")
  PLOT.DEFAULTS <- "workflow/scripts/plotDefaults.R"
  SPECIES <- "Bacteroides_finegoldii_57739"
  INDIR <- paste0("workflow/out/midasOutput/snps/HouseholdTransmission-Stool/",SPECIES,"/")
  OUTDIR <- paste0("workflow/report/calculateFixedDifferences/",SPECIES,"/")
  ifelse(!dir.exists(OUTDIR),dir.create(OUTDIR, recursive=TRUE),FALSE)
  HOUSEHOLD <- "XB"
} else{
  # List snakemake parameters required for running script.
  # snakemake@params[["plotDefaults"]] Path to plot default themes and palettes
  # snakemake@params[["species"]] Name of species to be analyzed
  # snakemake@params[["indir"]] Path to input directory containing MIDAS snp output files
  # snakemake@params[["outdir"]] Path to output directory for plots and small output files
  # snakemake@params[["household"]] Name of household to be analyzed
  PLOT.DEFAULTS <- snakemake@params[["plotDefaults"]]
  SPECIES <- snakemake@params[["species"]]
  INDIR <- paste0(snakemake@params[["indir"]],SPECIES,"/")
  OUTDIR <- paste0(snakemake@params[["outdir"]], SPECIES, "/")
  ifelse(!dir.exists(OUTDIR),dir.create(OUTDIR, recursive=TRUE),FALSE)
  HOUSEHOLD <- snakemake@params[["household"]]
}


# Import SNP and coverage data and combine into a single, filtered dataframe --------

dataFilterCovSamples <- importSNPsFilterSites(
  snpsDepthPath=paste0(INDIR,"snps_depth.txt.gz"),
  snpsFreqPath=paste0(INDIR,"snps_freq.txt.gz"),
  snpsInfoPath = paste0(INDIR,"snps_info.txt.gz"),
  outCoveragePath = paste0(OUTDIR,"sampleMedianCoverage.txt")
)
gc()


# Extract only the SNPs that belong to the household of interest. ---------

# Extract samples from the household of interest.
dataHousehold <- dataFilterCovSamples %>%
  filter(grepl(HOUSEHOLD, sample))
# Shorten the sample names to save space.
dataHousehold <- dataHousehold %>%
  mutate(sample=substr(sample, nchar(sample)-6, nchar(sample)))
# Export the resulting dataframe.
saveRDS(dataHousehold, paste0(INDIR, HOUSEHOLD, ".rds"))
