# For each species in each subject, annotate whether
# the species is shared with a cohabiting subject before abx,
# after abx, during the main study, and during follow-up sampling.
# Also annotate species sharing within the same subject
# before abx, after abx, during the main study, and during follow-up sampling.
library(tidyverse)
library(foreach)

# Import metadata about the study.
source("workflow/analysis/background/background.R")
# Set the output directory.
OUTDIR <- "workflow/analysis/annotateSpeciesTrajectoriesSharedStrains/out/"
LIMITOFDETECTION <- 1e-3


# Import the species abundance data.
source("workflow/analysis/generateSpeciesAbundances/loadSpeciesAbundancesFiltered.R")

# Summarize within-household species sharing.
# For each species in a subject that's present above the limit of detection,
# determine whether it is detected at all in any cohabiting partner.
# First, summarize the maximum pre-abx, post-abx, main study, and follow-up
# relative abundances for each species.
speciesAbundancesSummary <- speciesAbundances %>%
  group_by(hh, subject, species_id) %>%
  summarize(maxRelAbundance=max(relative_abundance),
            maxPreAbxRelAbundance=ifelse(sum(timepoint<30)>0,
                                         max(relative_abundance[timepoint<30]), 0),
            maxPostAbxRelAbundance=ifelse(sum(timepoint>34 & timepoint<75)>0,
                                          max(relative_abundance[timepoint>34 & timepoint<75]), 0),
            maxMainStudyRelAbundance=ifelse(sum(timepoint<75)>0,
                                            max(relative_abundance[timepoint<75]), 0),
            maxFollowupRelAbundance=ifelse(sum(timepoint>75)>0,
                                            max(relative_abundance[timepoint>75]), 0),
            finalRelAbundance=relative_abundance[timepoint==max(timepoint)])
# For each species in each subject present above the limit of detection, 
# determine if the species is detected
# in each period of the study in other cohabiting subjects.
speciesAbundancesHhComparison <- 
  foreach(isubject=subjectsX, .combine="rbind") %do% {
    # Extract the relative abundance summaries for species present
    # above the limit of detection in the main subject.
    focalSubjectSpecies <- speciesAbundancesSummary %>%
      filter(maxRelAbundance>LIMITOFDETECTION, subject==isubject)
    # Extract the relative abundance summaries for species present
    # in the cohabiting partners. When there are multiple cohabiting partners,
    # calculate the maximum relative abundance in each time period
    # within the household as a whole.
    cohabitingSubjectsSpecies <- speciesAbundancesSummary %>%
      filter(subject!=isubject, hh==substr(isubject,1,2)) %>%
      group_by(species_id) %>%
      summarize(maxRelAbundanceHh=max(maxRelAbundance),
                maxPreAbxRelAbundanceHh=max(maxPreAbxRelAbundance),
                maxPostAbxRelAbundanceHh=max(maxPostAbxRelAbundance),
                maxMainStudyRelAbundanceHh=max(maxMainStudyRelAbundance),
                maxFollowupRelAbundanceHh=max(maxFollowupRelAbundance))
    # Combine the relative abundance data from the focal subject
    # and cohabiting subjects.
    focalSubjectSpecies <- left_join(focalSubjectSpecies, cohabitingSubjectsSpecies,
                                     by=c("species_id"))
    # Annotate species sharing in the cohabiting subjects.
    focalSubjectSpeciesSharing <- focalSubjectSpecies %>% ungroup() %>%
      mutate(speciesSharedHhAnyTimepointPresence=(maxRelAbundanceHh>0),
             speciesSharedHhPreAbxPresence=(maxPreAbxRelAbundanceHh>0),
             speciesSharedHhPostAbxPresence=(maxPostAbxRelAbundanceHh>0),
             speciesSharedHhMainStudyPresence=(maxMainStudyRelAbundanceHh>0),
             speciesSharedHhFollowupPresence=(maxFollowupRelAbundanceHh>0),
             speciesSharedHhAnyTimepointLOD=(maxRelAbundanceHh>LIMITOFDETECTION),
             speciesSharedHhPreAbxLOD=(maxPreAbxRelAbundanceHh>LIMITOFDETECTION),
             speciesSharedHhPostAbxLOD=(maxPostAbxRelAbundanceHh>LIMITOFDETECTION),
             speciesSharedHhMainStudyLOD=(maxMainStudyRelAbundanceHh>LIMITOFDETECTION),
             speciesSharedHhFollowupLOD=(maxFollowupRelAbundanceHh>LIMITOFDETECTION)) %>%
      dplyr::select(subject, species_id, maxRelAbundance, starts_with("speciesShared"))
    return(focalSubjectSpeciesSharing)
  }
write.table(speciesAbundancesHhComparison, paste0(OUTDIR, "speciesSharing-hh.txt"),
            row.names=FALSE, quote=FALSE)

# Also compare species abundances within the same subject.
# For each species that exceeds the limit of detection at any timepoint,
# determine whether it is detected at all and whether it is detected
# above the limit of detection in each sampling period of interest.
speciesAbundancesSubjectSummary <- speciesAbundancesSummary %>%
  filter(maxRelAbundance>LIMITOFDETECTION) %>% ungroup() %>%
  mutate(speciesSubjectAnyTimepointPresence=(maxRelAbundance>0),
         speciesSubjectPreAbxPresence=(maxPreAbxRelAbundance>0),
         speciesSubjectPostAbxPresence=(maxPostAbxRelAbundance>0),
         speciesSubjectMainStudyPresence=(maxMainStudyRelAbundance>0),
         speciesSubjectFollowupPresence=(maxFollowupRelAbundance>0),
         speciesSubjectFinalSeqPresence=(finalRelAbundance>0),
         speciesSubjectAnyTimepointLOD=(maxRelAbundance>LIMITOFDETECTION),
         speciesSubjectPreAbxLOD=(maxPreAbxRelAbundance>LIMITOFDETECTION),
         speciesSubjectPostAbxLOD=(maxPostAbxRelAbundance>LIMITOFDETECTION),
         speciesSubjectMainStudyLOD=(maxMainStudyRelAbundance>LIMITOFDETECTION),
         speciesSubjectFollowupLOD=(maxFollowupRelAbundance>LIMITOFDETECTION),
         speciesSubjectFinalSeqLOD=(finalRelAbundance>LIMITOFDETECTION)) %>%
  dplyr::select(subject, species_id, maxRelAbundance, starts_with("speciesSubject"))
write.table(speciesAbundancesSubjectSummary, paste0(OUTDIR, "speciesSharing-subject.txt"),
            row.names=FALSE, quote=FALSE)
