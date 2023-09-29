# This script takes the MIDAS output files, filters out sites with abnormal coverage,
# and performs strain "fishing" by looking for the presence of private SNPs
# in all other samples.

library(data.table)
library(tidyverse)
library(cowplot)
library(foreach)
library(forcats)
library(R.utils)

# If TEST is true, then run on default files
# rather than snakemake inputs.
TEST <- FALSE
if(TEST){
  setwd("Z:/oak/stanford/groups/relman/users/kxue/household-transmission-mgx")
  PLOT.DEFAULTS <- "workflow/scripts/plotDefaults.R"
  SPECIES <- "Bacteroides_stercoris_56735"
  INDIR <- paste0("workflow/out/midasOutput/snps/HouseholdTransmission-Stool/",SPECIES,"/")
  OUTDIR <- paste0("workflow/report/performStrainFishing/",SPECIES,"/")
  ifelse(!dir.exists(OUTDIR),dir.create(OUTDIR, recursive=TRUE),FALSE)
  HMPSNPFREQS <- paste0("workflow/out/midasOutput/HMP/snp_prevalences/",SPECIES,".txt.gz")
} else{
  # List snakemake parameters required for running script.
  # snakemake@params[["plotDefaults"]] Path to plot default themes and palettes
  # snakemake@params[["species"]] Name of species to be analyzed
  # snakemake@params[["indir"]] Path to input directory containing MIDAS snp output files
  # snakemake@params[["outdir"]] Path to output directory for plots and small output files
  # snakemake@params[["HMPsnpsdir"]] Path to directory with files describing HMP SNP frequencies
  PLOT.DEFAULTS <- snakemake@params[["plotDefaults"]]
  SPECIES <- snakemake@params[["species"]]
  INDIR <- paste0(snakemake@params[["indir"]],SPECIES,"/")
  OUTDIR <- paste0(snakemake@params[["outdir"]], SPECIES, "/")
  ifelse(!dir.exists(OUTDIR),dir.create(OUTDIR, recursive=TRUE),FALSE)
  HMPSNPFREQS <- paste0(snakemake@params[["HMPsnpsdir"]],SPECIES,".txt.gz")
}

# Import helper functions for parsing the MIDAS output files.
source("workflow/scripts/helperFunctions.R")


# Import SNP and coverage data and combine into a single, filtered dataframe --------

# Import SNP data and filter out sites with abnormally low coverage.
dataFilterCovSamples <- importSNPsFilterSites(
  snpsDepthPath=paste0(INDIR,"snps_depth.txt.gz"),
  snpsFreqPath=paste0(INDIR,"snps_freq.txt.gz"),
  snpsInfoPath = paste0(INDIR,"snps_info.txt.gz"),
  outCoveragePath = paste0(OUTDIR,"sampleMedianCoverage.txt")
)
gc()

# Parse the sample names.
dataFilterCovSamples <- dataFilterCovSamples %>%
  mutate(sample=substr(sample,29,36),
         household=substr(sample,1,2),
         subject=substr(sample,1,3),
         timepoint=as.numeric(substr(sample,5,7)))

# Calculate the median coverage per sample.
dataMedianCoverage <- dataFilterCovSamples %>%
  group_by(sample) %>%
  summarize(cov=median(coverage)) %>%
  arrange(cov)

# Identify quasi-phaseable samples --------------------------------------------------

# Identify quasi-phaseable samples.
# Calculate the number of SNPs per sample that fall between 20%-80% frequency.
# When the number of SNPs in this range is on the order of 1% of the genome,
# the sample is likely to have complex strain structure.
# Use the dataframe in which low-coverage sites and samples have been filtered out.
dataIntermediateSNPs <- dataFilterCovSamples %>%
  filter(freq>0.2 & freq<0.8) %>%
  group_by(sample) %>% 
  summarize(numIntermediates=n())

# Classify samples that have few intermediate-frequency SNPs as quasi-phaseable.
# If the number of intermediate-frequency SNPs exceeds 0.1% of the genome length,
# then the sample is not quasi-phaseable.
samplesQP <- (dataIntermediateSNPs %>%
                filter(numIntermediates< (0.001*max(dataFilterCovSamples$site_id))))$sample


# Identify the global major allele at each site and polarize frequ --------

# Determine which allele is the global major and minor allele.
# First, choose one timepoint per subject.
# Choose the earliest timepoint that has coverage of at least 5,
# and exclude the subject if no timepoints have coverage of at least 5.
# Exclude subjects labeled as WWW (unknown).
oneSamplePerSubject <- (dataMedianCoverage %>%
                          mutate(subject=substr(sample,1,3), timepoint=as.numeric(substr(sample,5,7))) %>%
                          filter(cov>=5, subject !="WWW") %>% group_by(subject) %>% top_n(-1, timepoint))$sample
# Using one representative sample per subject, determine for each SNP site
# if the global major allele is 0 or 1.
# The global major allele is the allele that is the majority allele in a majority
# of the subjects in the dataset.
globalMajorAlleles <- dataFilterCovSamples %>%
  filter(sample %in% oneSamplePerSubject) %>%
  # Identify the major allele in each subject as 0 or 1.
  mutate(majorAllele=ifelse(freq<0.5,0,1)) %>%
  # Exclude sites with insufficient coverage.
  filter(!is.na(freq)) %>%
  # For each site, tally the number of subjects that have 1 as the major allele.
  group_by(site_id) %>%
  summarize(numSubjectsMajorAllele1=sum(majorAllele), numSubjectsMajorAllele=n(),
            freqMajorAllele1=numSubjectsMajorAllele1/numSubjectsMajorAllele) %>%
  # Require that at least 75% of subjects have adequate coverage at a site
  # before calling the global major allele.
  # Some parts of the genome may only be present in a small number of subjects.
  filter(numSubjectsMajorAllele < 0.75*length(oneSamplePerSubject)) %>%
  # Identify the global major allele, i.e. the allele that is the majority allele
  # in a majority of subjects.
  mutate(globalMajorAllele=ifelse(freqMajorAllele1<0.5,0,1))
# Extract vectors of sites at which the global major allele is 0 or 1.
globalMajorAllele0 <- (globalMajorAlleles %>% filter(globalMajorAllele==0))$site_id
globalMajorAllele1 <- (globalMajorAlleles %>% filter(globalMajorAllele==1))$site_id
# Create a globalFreq column in which allele frequencies are polarized such that
# 0 is the global majority allele and 1 is the global minority allele.
dataFilterCovSamples <- dataFilterCovSamples %>%
  mutate(globalFreq=ifelse(site_id %in% globalMajorAllele0, freq,
                    ifelse(site_id %in% globalMajorAllele1, 1-freq, NA)))


# Identify private SNP sites for this species -----------------------------

if(file.exists(HMPSNPFREQS)){
  
  # Import the list of SNP frequencies in the Human Microbiome Project for this species,
  # if available.
  HMPfrequencies <- fread(HMPSNPFREQS, header=TRUE, stringsAsFactors = FALSE, data.table=FALSE)
  
  # Import the site information for this species (a MIDAS output).
  dataInfo <- fread(paste0(INDIR,"snps_info.txt.gz"), header=TRUE, stringsAsFactors = FALSE, data.table=FALSE)
  
  # Bind site_id column of dataInfo to the HMPfrequencies dataset.
  # That is, link the site_id to the chromosome | location identification in the HMPfrequencies data.
  # Note that some sites in the HMPfrequencies dataset will have a site_id of N/A
  # if they were not called as variable sites in the household cohort.
  # However, this does not affect downstream analyses because these sites are
  # presumably invariant in my dataset.
  HMPfrequencies <- left_join(HMPfrequencies,
                              dataInfo %>% dplyr::select(site_id, ref_id, ref_pos) %>%
                                dplyr::rename(Chromosome=ref_id, Location=ref_pos),
                              by=c("Chromosome","Location"))
  
  # Define public SNPs as sites at which the AltFreq is >0,
  # that is, there are non-zero subjects in the HMP at which the global minor allele
  # is the major allele in that subject.
  publicSites <- sort(unique((HMPfrequencies %>% filter(AltFreq>0))$site_id))
  
  # Identify the private sites in the household cohort dataset.
  # These are sites that are not in the list of public sites from the HMP.
  privateSites <- sort(unique((dataFilterCovSamples %>% 
                                 filter(!(site_id %in% publicSites)))$site_id))
  
  # Iterate through each "bait" QP sample, identify private marker SNPs
  # (private sites in which the majority allele in the bait sample is the global minority allele),
  # and look for any evidence of variation at these private marker SNP sites
  # in the "query" samples.
  # Set a counter for QP samples.
  isamplesQP <- 0
  if(length(samplesQP)>0 & length(privateSites)>0){
    dataStrainFishing <- foreach(bait=samplesQP, .combine="rbind") %do% {
      # For a "bait" QP sample, identify private marker SNPs.
      # These private marker SNPs are sites that are 1) not in the list of public SNPs
      # and 2) have the global minor allele at high frequency in the bait.
      privateMarkerSNPs <- dataFilterCovSamples %>%
        filter(sample==bait, site_id %in% privateSites, globalFreq>0.8)
      
      # For all "query" samples, analyze the private marker SNP sites
      # in the "bait" sample that have non-zero frequency in the "query" sample.
      freqsPrivateMarkerSNPs <- dataFilterCovSamples %>%
        filter(site_id %in% privateMarkerSNPs$site_id,
               !is.na(globalFreq))
      # Determine how many private marker SNPs are detected in each "query" sample.
      numPrivateMarkerSNPs <- freqsPrivateMarkerSNPs %>%
        group_by(sample) %>%
        summarize(numSitesDetected=sum(globalFreq>0),
                  numSitesAvailable=n(),
                  meanSNPfreq=mean(globalFreq),
                  medianSNPfreq=median(globalFreq)) %>%
        mutate(bait=bait)
      isamplesQP <- isamplesQP+1
      print(paste0(isamplesQP," of ",length(samplesQP), " samples complete."))
      if(nrow(privateMarkerSNPs)>0){
        return(numPrivateMarkerSNPs)
      } else{
        return()
      }
    }
  }
  
  # Export the strain fishing data.
  write.table(dataStrainFishing, gzfile(paste0(OUTDIR,"strainFishing.txt.gz")),
              quote=FALSE, row.names = FALSE)
  
}

# Touch a file to verify that the script finished running.
file.create(paste0(OUTDIR,"done.txt"))