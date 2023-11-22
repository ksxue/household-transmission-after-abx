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
  setwd("Z:/oak/stanford/groups/relman/users/kxue/household-transmission-mgx")
  PLOT.DEFAULTS <- "workflow/scripts/plotDefaults.R"
  SPECIES <- "Bacteroides_stercoris_56735"
  INDIR <- paste0("workflow/out/midasOutput/snps/HouseholdTransmission-Stool/",SPECIES,"/")
  OUTDIR <- paste0("workflow/report/calculateFixedDifferences/",SPECIES,"/")
  ifelse(!dir.exists(OUTDIR),dir.create(OUTDIR, recursive=TRUE),FALSE)
} else{
  # List snakemake parameters required for running script.
  # snakemake@params[["plotDefaults"]] Path to plot default themes and palettes
  # snakemake@params[["species"]] Name of species to be analyzed
  # snakemake@params[["indir"]] Path to input directory containing MIDAS snp output files
  # snakemake@params[["outdir"]] Path to output directory for plots and small output files
  PLOT.DEFAULTS <- snakemake@params[["plotDefaults"]]
  SPECIES <- snakemake@params[["species"]]
  INDIR <- paste0(snakemake@params[["indir"]],SPECIES,"/")
  OUTDIR <- paste0(snakemake@params[["outdir"]], SPECIES, "/")
  ifelse(!dir.exists(OUTDIR),dir.create(OUTDIR, recursive=TRUE),FALSE)
}


# Import SNP and coverage data and combine into a single, filtered dataframe --------

dataFilterCovSamples <- importSNPsFilterSites(
  snpsDepthPath=paste0(INDIR,"snps_depth.txt.gz"),
  snpsFreqPath=paste0(INDIR,"snps_freq.txt.gz"),
  snpsInfoPath = paste0(INDIR,"snps_info.txt.gz"),
  outCoveragePath = paste0(OUTDIR,"sampleMedianCoverage.txt")
)
gc()

# Plot the site-frequency spectrum ----------------------------------------

# Set cowplot theme.
theme_set(theme_cowplot())
# Import plot defaults.
source(PLOT.DEFAULTS)

numSamples=n_distinct(dataFilterCovSamples$sample)
p <- dataFilterCovSamples %>%
  # Calculate minor allele frequency
  mutate(minorFreq=ifelse(freq>0.5,1-freq,freq)) %>%
  filter(minorFreq>0.01) %>%
  ggplot() +
  geom_histogram(aes(x=freq)) +
  xlim(0,0.5) +
  facet_wrap(~sample, scales="free_y", ncol=8) +
  guides(fill=FALSE) +
  ggtitle(SPECIES) +
  DEFAULTS.THEME_ALL
save_plot(paste0(OUTDIR,"SFS.png"), p, 
          ncol=2, nrow=max(0.25*ceiling(numSamples/8),1), limitsize=FALSE)

# Calculate the number of fixed differences between populations ----------------------------------------

# List the distinct samples.
samples <- unique(dataFilterCovSamples$sample)
if(length(samples)>1){
  # Calculate the number of fixed differences between populations,
  # also known as popANI as defined by Olm et al.
  dataFilterCovSamplesWide <- dataFilterCovSamples %>%
    dplyr::select(site_id, sample, freq) %>%
    spread(sample, freq, fill=NA)
  samplePairs <- combn(c(samples), 2)
  iSamplePairs <- 0
  pairwisePopANI <-
    foreach(sample1=samplePairs[1,],
            sample2=samplePairs[2,], .combine="rbind") %do%
    {
      distances <- dataFilterCovSamplesWide %>%
        mutate(diff=abs(!!as.name(sample1)-!!as.name(sample2)))
      iSamplePairs <- iSamplePairs+1
      if(iSamplePairs %% 1000 == 0){
        print(paste0(iSamplePairs," of ",length(iSamplePairs), " samples complete.")) 
      }
      return(c(sample1, sample2, 
               nrow(distances %>% filter(diff==1)),
               nrow(distances %>% filter(!is.na(diff)))))
    }
  # Rearrange the dataframe of fixed differences if there is only one sample pair.
  if(ncol(samplePairs)==1){
    fixedDiffs <- as.data.frame(t(pairwisePopANI))
  }else{
    fixedDiffs <- as.data.frame(pairwisePopANI)
  }
  colnames(fixedDiffs) <- c("sample1","sample2","fixedDiffs","sitesCompared")
  fixedDiffs$fixedDiffs <- as.integer(fixedDiffs$fixedDiffs)
  fixedDiffs$sitesCompared <- as.integer(fixedDiffs$sitesCompared)
  fixedDiffs$sample1 <- as.character(fixedDiffs$sample1)
  fixedDiffs$sample2 <- as.character(fixedDiffs$sample2)
  # Arrange the samples in each pair alphabetically.
  fixedDiffs <- fixedDiffs %>%
    mutate(orderedsample1=ifelse(as.character(sample1)<as.character(sample2),
                                 sample1,sample2),
           orderedsample2=ifelse(as.character(sample1)<as.character(sample2),
                                 sample2,sample1),
           sample1=orderedsample1, sample2=orderedsample2) %>%
    dplyr::select(-orderedsample1, -orderedsample2)
  # Export the list of fixed differences.
  write.table(fixedDiffs %>% arrange(sample1, sample2), 
              paste0(OUTDIR,"fixedDiffs_filterCoverage.txt"),
              quote=FALSE, row.names = FALSE) 
}


# Touch a file to verify that the script finished running.
file.create(paste0(OUTDIR,"done.txt"))
