library(tidyverse)

# Import metadata about the study.
source("workflow/analysis/background/background.R")
# Import plot defaults.
source("workflow/analysis/plotDefaults.R")
# Set the output directory.
OUTDIR <- "workflow/analysis/inferXBAstrainDynamicsPersistenceRecolonization/out/"

# Write a function that infers the SNPs that distinguish the initial strains of a species
# in subjects XBA and XBB and then extracts their median frequencies over time.
extractDistinguishingSNPs <- function(iSpecies){
  # Import the SNP data from household XB for the species specified.
  dataFreqs <- readRDS(paste0("workflow/out/midasOutput/snps/HouseholdTransmission-Stool/",
                                 iSpecies, "/XB.rds"))
  dataFreqs <- dataFreqs %>%
    mutate(subject=substr(sample,1,3), timepoint=as.integer(substr(sample,5,7)))
  # Calculate the median coverage for each sample.
  # Also calculate the number of SNPs with intermediate frequency,
  # i.e. between 0.2 and 0.8.
  # Define samples as QP (quasi-phaseable) if the number of intermediate SNPs
  # is less than 0.1% the genome length.
  dataCov <- dataFreqs %>%
    group_by(sample, subject, timepoint) %>%
    summarize(medianCov=median(coverage),
              numIntermediateSNPs=sum(!is.na(freq) & freq>0.2 & freq<0.8),
              genomeLength=max(site_id),
              QP=(numIntermediateSNPs<0.0015*genomeLength))
  # Plot the site-frequency spectrum of samples with coverage >5.
  p_SFS <- dataFreqs %>%
    filter(sample %in% (dataCov %>% filter(medianCov>5))$sample) %>%
    mutate(polarizedFreq=ifelse(freq<0.5,freq,1-freq)) %>%
    filter(polarizedFreq!=0) %>%
    ggplot() +
    geom_histogram(aes(x=polarizedFreq, fill=factor(subject)), binwidth=0.05) +
    facet_wrap(~sample, ncol=10, scales="free") +
    scale_fill_brewer(palette="Set1", name="Subject") +
    xlab("Minor allele frequency") + ylab("Number of sites") +
    DEFAULTS.THEME_PRINT
  save_plot(paste0(OUTDIR, iSpecies, "-SFS.png"), p_SFS, 
            base_width=10, base_height=n_distinct((dataCov %>% filter(medianCov>5))$sample)/10*0.8)
    
  # Calculate the number of subjects with QP samples with coverage >5
  # at at least one timepoint before antibiotics and after day 400.
  dataCovSummary <- dataCov %>%
    group_by(subject) %>%
    summarize(numPreAbx=sum(medianCov>5 & timepoint<=29 & QP), 
              numPostSwitch=sum(medianCov>5 & timepoint>400 & QP))
  nSubjects <- nrow(dataCovSummary %>% 
                      filter(numPreAbx>0, numPostSwitch>0))
  # Identify the QP samples in this dataset.
  QPsamples <- (dataCov %>% filter(QP))$sample
  # Continue only if both subjects have high-coverage samples before and after antibiotics.
  if(nSubjects==2){
    
    # Identify sites that have high max frequency (>0.8) in one subject 
    # and low max frequency (<0.2) in the other subject before antibiotics.
    # These sites would be quasi-phased as distinct sites.
    # Focus only on QP samples when identifying these sites.
    diffSitesSubjects <- dataFreqs %>%
      filter(timepoint<=29, !is.na(freq), sample %in% QPsamples) %>%
      group_by(subject, site_id) %>%
      summarize(maxFreq=max(freq)) %>%
      ungroup() %>% 
      pivot_wider(names_from=subject, values_from=maxFreq) %>%
      filter((XBA>0.8 & XBB<0.2) | (XBA<0.2 & XBB>0.8))
    # Extract the frequencies of alleles at distinguishing sites.
    diffSitesSubjectsFreqs <- dataFreqs %>%
      filter(site_id %in% diffSitesSubjects$site_id)
    # Polarize the site frequencies so that the XBA allele always has frequency 0.
    diffSitesSubjectsFreqs <- diffSitesSubjectsFreqs %>%
      mutate(polarizedFreq=
               ifelse(site_id %in% (diffSitesSubjects %>% filter(XBA>0.8))$site_id, 1-freq, freq))
    # Export the list of sites that are different between subjects.
    write.table(diffSitesSubjects, paste0(OUTDIR, iSpecies, "-diffSitesSubjects.txt"),
                row.names=FALSE, quote=FALSE)
    # Export the frequencies of these sites.
    write.table(diffSitesSubjectsFreqs, gzfile(paste0(OUTDIR, iSpecies, "-diffSitesSubjects-freqs.txt.gz")),
                row.names=FALSE, quote=FALSE)
    
    # Identify sites that have high mean frequency (>0.8) in the pre-abx XBA timepoints
    # and low mean frequency (<0.2) in the samples after the switch (or vice versa).
    diffSitesTimePeriodsXBA <- dataFreqs %>%
      filter(subject=="XBA", !is.na(freq),
             (timepoint<=29 | timepoint>400), sample %in% QPsamples) %>%
      mutate(timePeriod=ifelse(timepoint<=29, "before", "after")) %>%
      group_by(timePeriod, site_id) %>%
      summarize(maxFreq=max(freq)) %>%
      ungroup() %>% pivot_wider(names_from=timePeriod, values_from=maxFreq) %>%
      filter((before>0.8 & after<0.2) | (before<0.2 & after>0.8))
    # Extract the frequencies of alleles at sites that vary after the switch.
    diffSitesTimePeriodsXBAFreqs <- dataFreqs %>%
      filter(site_id %in% diffSitesTimePeriodsXBA$site_id)
    # Polarize the site frequencies so that the XBA allele always has frequency 0.
    diffSitesTimePeriodsXBAFreqs <- diffSitesTimePeriodsXBAFreqs %>%
      mutate(polarizedFreq=
               ifelse(site_id %in% (diffSitesTimePeriodsXBA %>% filter(before>0.8))$site_id, 1-freq, freq))
    # Export the list of sites that are different between subjects.
    write.table(diffSitesTimePeriodsXBA, paste0(OUTDIR, iSpecies, "-diffSitesTimePeriodsXBA.txt"),
                row.names=FALSE, quote=FALSE)
    # Export the frequencies of these sites.
    write.table(diffSitesTimePeriodsXBAFreqs, gzfile(paste0(OUTDIR, iSpecies, "-diffSitesTimePeriodsXBA-freqs.txt.gz")),
                row.names=FALSE, quote=FALSE)
    
    # Identify sites that consistently have intermediate frequencies
    # in all of the non-QP samples.
    # These are the sites that likely distinguish co-colonizing strains.
    diffSitesCocolonizingStrains <- dataFreqs %>%
      filter(!is.na(freq), !(sample %in% QPsamples), freq!=0, freq!=1) %>%
      group_by(subject, site_id) %>%
      summarize(numTimepointsIntermediateFreq=n_distinct(timepoint)) %>%
      ungroup() %>% group_by(subject) %>%
      filter(numTimepointsIntermediateFreq==max(numTimepointsIntermediateFreq)) %>%
      dplyr::select(subject, site_id)
    diffSitesCocolonizingStrainsXBA <- diffSitesCocolonizingStrains %>%
      filter(subject=="XBA")
    diffSitesCocolonizingStrainsXBB <- diffSitesCocolonizingStrains %>%
      filter(subject=="XBB")
    # Extract the frequencies of SNPs at sites that distinguish co-colonizing strains.
    diffSitesCocolonizingStrainsFreqs <- dataFreqs %>%
      filter((site_id %in% diffSitesCocolonizingStrainsXBA$site_id) | 
               (site_id %in% diffSitesCocolonizingStrainsXBB$site_id)) %>%
      mutate(polarizedFreq=ifelse(freq<0.5,freq,1-freq))
    # Export the list of sites that distinguish cocolonizing strains.
    write.table(diffSitesCocolonizingStrains, paste0(OUTDIR, iSpecies, "-diffSitesCocolonizingStrains.txt"),
                row.names=FALSE, quote=FALSE)
    # Export the frequencies of these sites.
    write.table(diffSitesCocolonizingStrainsFreqs, gzfile(paste0(OUTDIR, iSpecies, "-diffSitesCocolonizingStrains-freqs.txt.gz")),
                row.names=FALSE, quote=FALSE)
    
  }
    
}

postBfragilisSpecies <- c("Bacteroides_fragilis_54507",
                          "Eggerthella_lenta_56463",
                          "Bifidobacterium_longum_57796",      
                          "Bacteroides_vulgatus_57955",
                          "Bacteroides_faecis_58503",
                          "Bacteroides_cellulosilyticus_58046",
                          "Bacteroides_finegoldii_57739",
                          "Bacteroides_ovatus_58035",
                          "Bacteroides_uniformis_57318",       
                          "Bacteroides_xylanisolvens_57185",
                          "Bacteroides_intestinalis_61596")
sapply(postBfragilisSpecies, extractDistinguishingSNPs)
