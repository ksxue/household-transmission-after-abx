library(tidyverse)


# Import metadata about the study.
source("workflow/analysis/background/background.R")
# Import plot defaults.
source("workflow/analysis/plotDefaults.R")
# Set the output directory.
OUTDIR <- "workflow/analysis/inferXBAstrainDynamicsPersistenceRecolonization/out/"

# Infer the overall species abundance dynamics for XBA and XBB.
source("workflow/analysis/generateSpeciesAbundances/loadSpeciesAbundancesFiltered.R")

# List the species of interest that undergo major changes in relative abundance.
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

speciesAbundancesXB <- speciesAbundances %>%
  filter(species_id %in% postBfragilisSpecies,
         subject %in% c("XBA", "XBB"))

# Infer the strain dynamics for B. uniformis ------------------------------

# Infer the species dynamics for B. uniformis.
speciesAbundancesXBBuniformis <- speciesAbundancesXB %>%
  filter(species_id=="Bacteroides_uniformis_57318") %>%
  dplyr::select(sample, species_id, subject, timepoint, relative_abundance)

# Import SNP files for B. uniformis for sites that distinguish XBA and XBB
dataXBBuniformis <- 
  read.table("workflow/analysis/inferXBAstrainDynamicsPersistenceRecolonization/out/Bacteroides_uniformis_57318-diffSitesSubjects-freqs.txt.gz",
             header=TRUE, stringsAsFactors = FALSE)
# Summarize the SNP dynamics for B. uniformis.
# Note that B. uniformis does not have appreciable strain structure (i.e. co-colonizing strains),
# and there is also not consistent evidence for SNP sweeps during this period.
dataXBBuniformisSummary <- dataXBBuniformis %>%
  filter(!is.na(freq)) %>%
  group_by(subject, timepoint) %>%
  summarize(strainFreq=median(polarizedFreq))
# Convert the average SNP dynamics for B. uniformis into strain dynamics.
strainAbundancesXBBuniformis <- dataXBBuniformisSummary %>%
  mutate(strain1=1-strainFreq, strain2=strainFreq) %>%
  dplyr::select(-strainFreq) %>%
  pivot_longer(c("strain1", "strain2"), names_to="strain", values_to="strainFreq")

# Combine the species and strain dynamics.
speciesStrainAbundancesXBBuniformis <-
  left_join(speciesAbundancesXBBuniformis,
            strainAbundancesXBBuniformis,
            by=c("subject", "timepoint")) %>%
  mutate(strainRelativeAbundance=strainFreq*relative_abundance) 
# At timepoints where the sequencing depth is not sufficient to infer strains,
# label the strain as "strainUnknown" to avoid making erroneous interpolations.
speciesStrainAbundancesXBBuniformis <- speciesStrainAbundancesXBBuniformis %>%
  mutate(strain=ifelse(is.na(strainRelativeAbundance), "strainUnknown", strain),
         strainRelativeAbundance=
           ifelse(is.na(strainRelativeAbundance), relative_abundance, strainRelativeAbundance)) %>%
  dplyr::select(subject, timepoint, species_id, strain, relative_abundance, strainFreq, strainRelativeAbundance) %>%
  group_by(subject) %>%
  complete(timepoint, strain, 
           fill=list(strainRelativeAbundance=0, strainFreq=0, species_id="Bacteroides_uniformis_57318")) %>%
  group_by(subject, timepoint) %>%
  mutate(relative_abundance=sum(strainRelativeAbundance))
# Export strain dynamics.
write.table(speciesStrainAbundancesXBBuniformis, 
            paste0(OUTDIR, "strainAbundances-Bacteroides_uniformis_57318.txt"),
            row.names=FALSE, quote=FALSE)


# Infer the strain dynamics for Bifidobacterium longum --------------------


# Infer the species dynamics for B. longum
speciesAbundancesXBBlongum <- speciesAbundancesXB %>%
  filter(species_id=="Bifidobacterium_longum_57796") %>%
  dplyr::select(sample, species_id, subject, timepoint, relative_abundance)

# Import SNP files for B. uniformis for sites that distinguish XBA and XBB
dataXBBlongum <- 
  read.table("workflow/analysis/inferXBAstrainDynamicsPersistenceRecolonization/out/Bifidobacterium_longum_57796-diffSitesSubjects-freqs.txt.gz",
             header=TRUE, stringsAsFactors = FALSE)
dataXBBlongumXBAXBBDistinguishingSites <- sort(unique(dataXBBlongum$site_id))

# Perform a further filtering step to remove sites that have highly variable frequencies,
# which are likely to result from mismapping and other errors.
dataXBBlongumHighQualitySNPs <- 
  unique((dataXBBlongum %>%
            group_by(site_id) %>%
            mutate(initFreqXBA=polarizedFreq[subject=="XBA" & timepoint==1],
                   initFreqXBB=polarizedFreq[subject=="XBB" & timepoint==1]) %>%
            filter(initFreqXBA==0, initFreqXBB==1))$site_id)
# Import the list of SNPs that may distinguish co-colonizing strains.
# Based on my review of these SNPs, it's not completely clear that they constitute a separate strain,
# but in any case, these sites that frequently have intermediate frequencies
# may present some issues for inferring strain dynamics,
# so I will remove them from the list of SNPs used to infer strain frequencies.
dataXBBlongumCocolonizing <- 
  read.table("workflow/analysis/inferXBAstrainDynamicsPersistenceRecolonization/out/Bifidobacterium_longum_57796-diffSitesCocolonizingStrains-freqs.txt.gz",
             header=TRUE, stringsAsFactors = FALSE)
dataXBBlongumCocolonizingSites <- sort(unique(dataXBBlongumCocolonizing$site_id))

# Count the number of sites that have a median frequency of 0 in XBA
# before antibiotics and a median frequency of 1 in XBA after the strain return.
# These sites also need to have median frequency of 1 in XBB as well.
# This is the number of sites that distinguishes the two strains.
dataXBBlongumXBAsweepingSites <- (dataXBBlongum %>%
                                    filter(!is.na(freq)) %>%
                                    group_by(site_id) %>%
                                    summarize(preAbxFreqXBA=median(polarizedFreq[timepoint<=29 & subject=="XBA"]),
                                              postAbxFreqXBA=median(polarizedFreq[timepoint>=400 & subject=="XBA"]),
                                              freqXBB=median(polarizedFreq[subject=="XBB"])) %>%
                                    filter(preAbxFreqXBA==0, postAbxFreqXBA==1, freqXBB==1))$site_id

length((dataXBBlongum %>%
          filter(!is.na(freq)) %>%
          group_by(site_id) %>%
          summarize(preAbxFreqXBA=median(polarizedFreq[timepoint<=29 & subject=="XBA"]),
                    postAbxFreqXBA=median(polarizedFreq[timepoint>=400 & subject=="XBA"]),
                    freqXBB=median(polarizedFreq[subject=="XBB"])) %>%
          filter(preAbxFreqXBA==0, freqXBB==1))$site_id)

# Summarize the SNP dynamics for B. longum between subjects.
# Note that B. longum may be comprised of multiple strains in XBA.
# There is also not consistent evidence for SNP sweeps during this period,
# since it also appears that there is a strain turnover.
dataXBBlongumSummary <- dataXBBlongum %>%
  filter(!is.na(freq), site_id %in% dataXBBlongumXBAsweepingSites) %>%
  group_by(subject, timepoint) %>%
  summarize(strainFreq=median(polarizedFreq),
            strainFreqLower=quantile(polarizedFreq,0.1),
            strainFreqUpper=quantile(polarizedFreq,0.9))
# Convert the average SNP dynamics for B. longum between hosts into strain dynamics.
strainAbundancesXBBlongum <- dataXBBlongumSummary %>%
  mutate(strain1=1-strainFreq, strain2=strainFreq) %>%
  dplyr::select(-strainFreq) %>%
  pivot_longer(c("strain1", "strain2"), names_to="strain", values_to="strainFreq")

# Combine the species and strain dynamics.
speciesStrainAbundancesXBBlongum <-
  left_join(speciesAbundancesXBBlongum,
            strainAbundancesXBBlongum,
            by=c("subject", "timepoint")) %>%
  mutate(strainRelativeAbundance=strainFreq*relative_abundance) 
# At timepoints where the sequencing depth is not sufficient to infer strains,
# label the strain as "strainUnknown" to avoid making erroneous interpolations.
speciesStrainAbundancesXBBlongum <- speciesStrainAbundancesXBBlongum %>%
  mutate(strain=ifelse(is.na(strainRelativeAbundance), "strainUnknown", strain),
         strainRelativeAbundance=
           ifelse(is.na(strainRelativeAbundance), relative_abundance, strainRelativeAbundance)) %>%
  dplyr::select(subject, timepoint, species_id, strain, relative_abundance, strainFreq, strainRelativeAbundance) %>%
  group_by(subject) %>%
  complete(timepoint, strain, 
           fill=list(strainRelativeAbundance=0, strainFreq=0, species_id="Bifidobacterium_longum_57796")) %>%
  group_by(subject, timepoint) %>%
  mutate(relative_abundance=sum(strainRelativeAbundance))
# Export strain dynamics.
write.table(speciesStrainAbundancesXBBlongum, 
            paste0(OUTDIR, "strainAbundances-Bifidobacterium_longum_57796.txt"),
            row.names=FALSE, quote=FALSE)


# Infer the strain dynamics for Bacteroides vulgatus --------------------

# Infer the species dynamics for B. vulgatus
speciesAbundancesXBBvulgatus <- speciesAbundancesXB %>%
  filter(species_id=="Bacteroides_vulgatus_57955") %>%
  dplyr::select(sample, species_id, subject, timepoint, relative_abundance)

# Import SNP files for B. vulgatus for sites that distinguish XBA and XBB
dataXBBvulgatus <- 
  read.table("workflow/analysis/inferXBAstrainDynamicsPersistenceRecolonization/out/Bacteroides_vulgatus_57955-diffSitesSubjects-freqs.txt.gz",
             header=TRUE, stringsAsFactors = FALSE)
dataXBBvulgatusXBAXBBDistinguishingSites <- sort(unique(dataXBBvulgatus$site_id))

# Summarize the SNP dynamics for B. vulgatus between subjects.
# Note that B. vulgatus may be comprised of multiple strains in XBB.
dataXBBvulgatusSummary <- dataXBBvulgatus %>%
  filter(!is.na(freq)) %>%
  group_by(subject, timepoint) %>%
  summarize(strainFreq=median(polarizedFreq))
# Convert the average SNP dynamics for B. vulgatus between hosts into strain dynamics.
strainAbundancesXBBvulgatus <- dataXBBvulgatusSummary %>%
  mutate(strain1=1-strainFreq, strain2=strainFreq) %>%
  dplyr::select(-strainFreq)

# Import the list of SNPs that may distinguish co-colonizing strains.
# XBB appears to consist of two co-colonizing strains.
dataXBBvulgatusCocolonizing <- 
  read.table("workflow/analysis/inferXBAstrainDynamicsPersistenceRecolonization/out/Bacteroides_vulgatus_57955-diffSitesCocolonizingStrains-freqs.txt.gz",
             header=TRUE, stringsAsFactors = FALSE)
dataXBBvulgatusCocolonizingSites <- sort(unique(dataXBBvulgatusCocolonizing$site_id))
# Polarize the sites that distinguish the co-colonizing strains in XBB
# so that the majority allele in XBA has frequency 1.
dataXBBvulgatusCocolonizingSitesXBAfreq <- dataXBBvulgatusCocolonizing %>%
  filter(subject=="XBA", !is.na(freq)) %>%
  group_by(site_id) %>% summarize(freq=median(freq))
dataXBBvulgatusCocolonizingSitesXBAfreq0 <- (dataXBBvulgatusCocolonizingSitesXBAfreq %>%
  filter(freq==0))$site_id
dataXBBvulgatusCocolonizingSitesXBAfreq1 <- (dataXBBvulgatusCocolonizingSitesXBAfreq %>%
                                               filter(freq==1))$site_id
dataXBBvulgatusCocolonizing <- dataXBBvulgatusCocolonizing %>%
  filter(!is.na(freq)) %>%
  mutate(polarizedFreqXBA=
           ifelse(site_id %in% dataXBBvulgatusCocolonizingSitesXBAfreq0, 1-freq,
           ifelse(site_id %in% dataXBBvulgatusCocolonizingSitesXBAfreq1, freq, NA)))
# Summarize the SNP dynamics for B. vulgatus co-colonizing strains in XBB.
# Note that B. vulgatus may be comprised of multiple strains in XBB.
dataXBBvulgatusCocolonizingSummary <- dataXBBvulgatusCocolonizing %>%
  filter(!is.na(polarizedFreqXBA)) %>%
  group_by(subject, timepoint) %>%
  summarize(strain2Freq=median(polarizedFreqXBA),
            strain2Freq25=quantile(polarizedFreqXBA,0.25),
            strain2Freq75=quantile(polarizedFreqXBA,0.75)) %>%
  mutate(strain3=1-strain2Freq)
# Add the frequencies of the co-colonizing strains in XBB onto the SNP dynamics for B. vulgatus.
# The SNPs whose frequencies are calculated here distinguish strain 2 from strain 3.
strainAbundancesXBBvulgatus <- 
  left_join(strainAbundancesXBBvulgatus, 
            dataXBBvulgatusCocolonizingSummary %>% dplyr::select(subject, timepoint, strain3),
            by=c("subject","timepoint")) %>% ungroup() %>%
  mutate(strain3=ifelse(is.na(strain3),0,strain3),
         strain2=(strain2-strain3))

# Import the list of SNPs that may represent sweeps in XBA.
dataXBAsweeps <- 
  read.table("workflow/analysis/inferXBAstrainDynamicsPersistenceRecolonization/out/Bacteroides_vulgatus_57955-diffSitesTimePeriodsXBA-freqs.txt.gz",
             header=TRUE, stringsAsFactors = FALSE)
# Retain only sites that have fewer than 3 NAs between both subjects.
dataXBAsweeps <- dataXBAsweeps %>%
  group_by(site_id) %>%
  filter(sum(is.na(polarizedFreq))<3)
# Summarize the mean frequency of the SNPs that sweep, which in this case appear to
# consist of a single haplotype.
dataXBAsweepsSummary <- dataXBAsweeps %>%
  filter(!is.na(polarizedFreq)) %>%
  group_by(subject, timepoint) %>%
  summarize(strain4=median(polarizedFreq))
# Add the frequencies of the SNP sweeps onto the strain dynamics for B. vulgatus.
# The SNPs whose frequencies are calculated here distinguish strain 4 from strain 1.
strainAbundancesXBBvulgatus <- 
  left_join(strainAbundancesXBBvulgatus, dataXBAsweepsSummary,
            by=c("subject","timepoint")) %>% ungroup() %>%
  mutate(strain4=ifelse(is.na(strain4),0,strain4),
         strain1=(strain1-strain4)) %>%
  mutate(strain1=ifelse(strain1<0,0,strain1))

# Tidy the strain abundances and normalize so that they sum to 1.
# Note that they may not initially sum to 1 because of variation in the measurement
# of the SNPs used to distinguish the strains.
strainAbundancesXBBvulgatus <- strainAbundancesXBBvulgatus %>%
  pivot_longer(c("strain1", "strain2","strain3","strain4"), 
               names_to="strain", values_to="strainFreq") %>%
  group_by(subject, timepoint) %>%
  mutate(strainFreq=strainFreq/sum(strainFreq))

# Combine the species and strain dynamics.
speciesStrainAbundancesXBBvulgatus <-
  left_join(speciesAbundancesXBBvulgatus,
            strainAbundancesXBBvulgatus,
            by=c("subject", "timepoint")) %>%
  mutate(strainRelativeAbundance=strainFreq*relative_abundance) 
# At timepoints where the sequencing depth is not sufficient to infer strains,
# label the strain as "strainUnknown" to avoid making erroneous interpolations.
speciesStrainAbundancesXBBvulgatus <- speciesStrainAbundancesXBBvulgatus %>%
  mutate(strain=ifelse(is.na(strainRelativeAbundance), "strainUnknown", strain),
         strainRelativeAbundance=
           ifelse(is.na(strainRelativeAbundance), relative_abundance, strainRelativeAbundance)) %>%
  dplyr::select(subject, timepoint, species_id, strain, relative_abundance, strainFreq, strainRelativeAbundance) %>%
  group_by(subject) %>%
  complete(timepoint, strain, 
           fill=list(strainRelativeAbundance=0, strainFreq=0, species_id="Bacteroides_vulgatus_57955")) %>%
  group_by(subject, timepoint) %>%
  mutate(relative_abundance=sum(strainRelativeAbundance))
# Export strain dynamics.
write.table(speciesStrainAbundancesXBBvulgatus, 
            paste0(OUTDIR, "strainAbundances-Bacteroides_vulgatus_57955.txt"),
            row.names=FALSE, quote=FALSE)


# Infer the strain dynamics for B faecis ----------------------------------


# Infer the species dynamics for B. faecis
speciesAbundancesXBBfaecis <- speciesAbundancesXB %>%
  filter(species_id=="Bacteroides_faecis_58503") %>%
  dplyr::select(sample, species_id, subject, timepoint, relative_abundance)

# Import SNP files for B. faecis for sites that distinguish XBA and XBB
dataXBBfaecis <- 
  read.table("workflow/analysis/inferXBAstrainDynamicsPersistenceRecolonization/out/Bacteroides_faecis_58503-diffSitesSubjects-freqs.txt.gz",
             header=TRUE, stringsAsFactors = FALSE)
# Retain only sites that have fewer than 20 NAs,
dataXBBfaecis <- dataXBBfaecis %>%
  group_by(site_id) %>%
  filter(sum(is.na(polarizedFreq))<20)
dataXBBfaecisXBAXBBDistinguishingSites <- sort(unique(dataXBBfaecis$site_id))

# Summarize the SNP dynamics for B. faecis between subjects.
# Group the SNPs into two sets based on whether or not their initial frequency was 0 in XBA.
dataXBBfaecisSummary <- dataXBBfaecis %>%
  filter(!is.na(freq)) %>%
  group_by(site_id) %>%
  mutate(SNPgroup=ifelse(polarizedFreq[timepoint==min(timepoint) & subject=="XBA"]==0,"SNPgroup1","SNPgroup2")) %>%
  ungroup() %>% group_by(subject, timepoint, SNPgroup) %>%
  summarize(strainFreq=median(polarizedFreq)) %>%
  ungroup() %>%
  pivot_wider(names_from=SNPgroup, values_from=strainFreq, values_fill=0)
# Convert the average SNP dynamics for B. faecis between hosts into strain dynamics.
strainAbundancesXBBfaecis <- dataXBBfaecisSummary %>%
  mutate(strain1=SNPgroup1, strain2=SNPgroup2-SNPgroup1, strain3=1-SNPgroup2) %>%
  dplyr::select(-SNPgroup1, -SNPgroup2)

# Import the list of SNPs that may represent sweeps in XBA.
dataXBAsweepsBfaecis <- 
  read.table("workflow/analysis/inferXBAstrainDynamicsPersistenceRecolonization/out/Bacteroides_faecis_58503-diffSitesTimePeriodsXBA-freqs.txt.gz",
             header=TRUE, stringsAsFactors = FALSE)
# After manually reviewing these SNPs, none appear to be likely candidates for sweeps.
# The SNPs that are more convincing are already included in SNP group 2 above.

# Tidy the strain abundances and normalize so that they sum to 1.
# Note that they may not initially sum to 1 because of variation in the measurement
# of the SNPs used to distinguish the strains.
strainAbundancesXBBfaecis <- strainAbundancesXBBfaecis %>%
  pivot_longer(c("strain1", "strain2","strain3"), 
               names_to="strain", values_to="strainFreq") %>%
  group_by(subject, timepoint) %>%
  mutate(strainFreq=strainFreq/sum(strainFreq))

# Combine the species and strain dynamics.
speciesStrainAbundancesXBBfaecis <-
  left_join(speciesAbundancesXBBfaecis,
            strainAbundancesXBBfaecis,
            by=c("subject", "timepoint")) %>%
  mutate(strainRelativeAbundance=strainFreq*relative_abundance) 
# At timepoints where the sequencing depth is not sufficient to infer strains,
# label the strain as "strainUnknown" to avoid making erroneous interpolations.
speciesStrainAbundancesXBBfaecis <- speciesStrainAbundancesXBBfaecis %>%
  mutate(strain=ifelse(is.na(strainRelativeAbundance), "strainUnknown", strain),
         strainRelativeAbundance=
           ifelse(is.na(strainRelativeAbundance), relative_abundance, strainRelativeAbundance)) %>%
  dplyr::select(subject, timepoint, species_id, strain, relative_abundance, strainFreq, strainRelativeAbundance) %>%
  group_by(subject) %>%
  complete(timepoint, strain, 
           fill=list(strainRelativeAbundance=0, strainFreq=0, species_id="Bacteroides_faecis_58503")) %>%
  group_by(subject, timepoint) %>%
  mutate(relative_abundance=sum(strainRelativeAbundance))
# Export strain dynamics.
write.table(speciesStrainAbundancesXBBfaecis, 
            paste0(OUTDIR, "strainAbundances-Bacteroides_faecis_58503.txt"),
            row.names=FALSE, quote=FALSE)

