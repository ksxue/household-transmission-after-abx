# These helper functions are meant for post-processing of MIDAS output.

library(data.table)
library(tidyverse)

# This function takes the merged MIDAS output for a given species,
# takes coverage and site info into account, generates a single dataframe,
# and filters out sites with abnormal coverage.
# It also exports the coverage per sample.
# Requires file paths for the MIDAS snps_depth, snps_freq, and snps_info files, used as input,
# as well as for an output path for the sample median coverage.
# Import SNP and coverage data and combine into a single, filtered dataframe.
# Repolarize SNP frequencies so that 0 is the reference allele and 1 is the alternate allele.
# Return the cleaned dataframe.
importSNPsFilterSites <- function(snpsDepthPath, snpsFreqPath, snpsInfoPath, outCoveragePath){
  
  # Import species relative abundances.
  dataCovRaw <- fread(snpsDepthPath, header=TRUE, stringsAsFactors = FALSE, data.table=FALSE)
  dataFreqRaw <- fread(snpsFreqPath, header=TRUE, stringsAsFactors = FALSE, data.table=FALSE)
  dataInfo <- fread(snpsInfoPath, header=TRUE, stringsAsFactors = FALSE, data.table=FALSE)
  
  # Tidy the data.
  dataCov <- gather(dataCovRaw, key="sample", value="coverage", -site_id)
  dataFreq <- gather(dataFreqRaw, key="sample", value="freq", -site_id)
  
  # Combine the coverage and frequency datasets.
  data <- left_join(dataCov, dataFreq, by=c("site_id","sample"))
  
  # Free memory from dataframes that are no longer needed.
  rm(dataCov)
  rm(dataFreq)
  
  # Calculate the median coverage per sample.
  dataMedianCoverage <- data %>%
    group_by(sample) %>%
    summarize(cov=median(coverage)) %>%
    arrange(cov)
  # Calculate the median coverage per sample after dropping sites with zero coverage.
  dataMedianCoverageNoZeroes <- data %>%
    filter(coverage>0) %>%
    group_by(sample) %>%
    summarize(covnozeroes=median(coverage)) %>%
    arrange(covnozeroes)
  # Merge the two metrics of coverage.
  dataMedianCoverage <- left_join(dataMedianCoverage,
                                  dataMedianCoverageNoZeroes, by="sample")
  # Export the coverage calculations for each sample to the desired location.
  write.table(dataMedianCoverage, outCoveragePath,
              quote=FALSE, row.names = FALSE)
  
  # Repolarize the SNPs so that 0 is the reference allele, and 1 is the alt allele,
  # based on the reference genome for the species.
  repolarize <- (dataInfo %>%
                   filter(ref_allele!=major_allele))$site_id
  data <- data %>%
    mutate(freq=ifelse(site_id %in% repolarize, 1-freq, freq))
  
  # Free memory from dataframes that are no longer needed.
  rm(dataInfo)
  
  # Exclude sites in each sample where the coverage is more than twofold
  # off from the genome-wide median for that sample.
  # If coverage does not fulfill this criteria, change the freq to NA.
  dataFilterCov <- data %>%
    group_by(sample) %>%
    mutate(freq=ifelse((coverage>median(coverage)/2)&(coverage<median(coverage)*2), 
                       freq, NA)) %>% ungroup()
  
  return(dataFilterCov)
}

