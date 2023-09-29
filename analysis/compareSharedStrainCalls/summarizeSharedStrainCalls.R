library(tidyverse)

# Import background metadata on the study.
source("workflow/analysis/background/background.R")
# Set the output directory.
OUTDIR <- "workflow/analysis/compareSharedStrainCalls/out/"


# Summarize strain sharing calls for samples from different subjec --------

# Import the list of consolidated strain sharing annotations.
sharedStrains <- read.table("workflow/analysis/compareSharedStrainCalls/out/strainSharingCalls-consolidated.txt",
                            header=TRUE, stringsAsFactors = FALSE) %>%
  dplyr::rename(species_id=species)

# Fill out the list of shared strain annotations so that each sample is
# annotated as sample 1 and sample 2.
sharedStrainsFull <- rbind(
  sharedStrains %>% dplyr::select(species_id, sample1, sample2, consensusNA),
  sharedStrains %>% dplyr::select(species_id, sample1, sample2, consensusNA) %>%
    dplyr::rename(sampletemp=sample1, sample1=sample2) %>% 
    dplyr::rename(sample2=sampletemp)
)
# Annotate the timepoints for the shared strain annotations.
sharedStrainsFull <- sharedStrainsFull %>%
  mutate(timepoint1=as.integer(substr(sample1,5,7)), 
         timepoint2=as.integer(substr(sample2,5,7)),
         subject1=substr(sample1,1,3),
         subject2=substr(sample2,1,3))

# For each species in each subject at each timepoint,
# determine if that species is ever shared with someone else in the household.
sharedStrainsByTimepointAllTimepoints <- sharedStrainsFull %>%
  filter(!is.na(consensusNA)) %>%
  group_by(species_id, sample1) %>%
  summarize(strainSharedAnyTimepoint=(sum(consensusNA)>0))
write.table(sharedStrainsByTimepointAllTimepoints,
            paste0(OUTDIR, "sameHousehold-sharedStrains-byTimepoint.txt"),
            row.names=FALSE, quote=FALSE)

# For each species in each subject across all timepoints,
# determine if that species is ever shared with someone else in the household
# and whether it is shared in different study periods (pre-abx, main, follow-up).
sharedStrainsBySpecies <- sharedStrainsFull %>%
  filter(!is.na(consensusNA)) %>%
  group_by(species_id, subject1) %>%
  summarize(strainSharedAnyTimepoint=(sum(consensusNA)>0),
         strainSharedPreAbx=ifelse(sum(timepoint1<30 & timepoint2<30)>0, 
                                   (sum(consensusNA[timepoint1<30 & timepoint2<30])>0), NA),
         strainSharedPostAbx=ifelse(sum(timepoint1>34 & timepoint1<75 & timepoint2>34 & timepoint2<75)>0, 
                                   (sum(consensusNA[timepoint1>34 & timepoint1<75 & timepoint2>34 & timepoint2<75])>0), NA),
         strainSharedMain=ifelse(sum(timepoint1<75 & timepoint2<75)>0, 
                                 (sum(consensusNA[timepoint1<75 & timepoint2<75])>0), NA),
         strainSharedFollowup=ifelse(sum(timepoint1>75 & timepoint2>75)>0, 
                                     (sum(consensusNA[timepoint1>75 & timepoint2>75])>0), NA),
         strainSharedPostAbxFollowup=ifelse(sum(timepoint1>34 & timepoint2>34)>0, 
                                            (sum(consensusNA[timepoint1>34 & timepoint2>34])>0), NA))
write.table(sharedStrainsBySpecies, paste0(OUTDIR, "sameHousehold-sharedStrains-subjectSummary.txt"),
            row.names=FALSE, quote=FALSE)



# Summarize strain sharing for samples from the same subject --------------

# Import the list of consolidated strain sharing annotations.
sharedStrainsSameSubject <- read.table("workflow/analysis/compareSharedStrainCalls/out/strainSharingCalls-consolidated-sameSubject.txt",
                            header=TRUE, stringsAsFactors = FALSE) %>%
  dplyr::rename(species_id=species)

# For each species in each subject at each timepoint, determine if the strain
# at that timepoint is shared with the initial or final strains from that subject.
sharedStrainsSameSubjectByTimepoint <- sharedStrainsSameSubject %>%
  filter(!is.na(consensusNA)) %>%
  filter(sample1 %in% samplesInitial | sample2 %in% samplesLastSequencedX |
           sample1 %in% samplesFinal | sample2 %in% samplesFinal) %>%
  mutate(focalSample=ifelse(sample1 %in% c(samplesInitial, samplesFinal), sample1,sample2),
         sample=ifelse(sample1==focalSample, sample2, sample1)) %>%
  dplyr::select(species_id, sample1, sample2, focalSample, sample, consensusNA) %>%
  arrange(focalSample, species_id, sample, consensusNA) %>%
  mutate(focalTimepoint=ifelse(focalSample %in% samplesInitial, "initial",
                        ifelse(focalSample %in% samplesFinal, "final", "lastSequenced")))
sharedStrainsSameSubjectByTimepointWide <- sharedStrainsSameSubjectByTimepoint %>%
  dplyr::select(-focalSample, -sample1, -sample2) %>%
  pivot_wider(names_from=focalTimepoint, values_from=consensusNA) %>%
  dplyr::rename(subjectSharedTimepointInitial=initial,
                subjectSharedTimepointFinal=final,
                subjectSharedTimepointLastSequenced=lastSequenced)
# Identify candidate strain turnovers among populations.
# Define a strain turnover as a switch from a strain being shared with the initial timepoint
# to a strain not being shared with the initial timepoint, with no switch back.
# If a strain turnover occurs, then annotate the first timepoint at which the new strain is detected.
sharedStrainsSameSubjectByTimepointWide <- sharedStrainsSameSubjectByTimepointWide %>%
  mutate(subject=substr(sample,1,3), timepoint=as.integer(substr(sample,5,7))) %>%
  group_by(subject, species_id) %>%
  filter(!is.na(subjectSharedTimepointInitial)) %>%
  mutate(strainSwitch=ifelse(timepoint==min(timepoint), !subjectSharedTimepointInitial,
                             subjectSharedTimepointInitial!=lag(subjectSharedTimepointInitial)),
         strainTurnover=(sum(strainSwitch)==1),
         timeOfStrainTurnover=ifelse(!strainTurnover, NA, min(timepoint[!subjectSharedTimepointInitial]))) %>%
  ungroup() %>% dplyr::select(-subject, -timepoint, -strainSwitch)

# Check the list of candidate strain turnovers (separate pipeline) and 
# import a list of verified strain turnovers based on the SNP trajectories.
verifiedStrainTurnovers <- read.table("workflow/analysis/scratch/230724-checkPreexistingStrainCalls/out/strainTurnovers-verified.txt",
                                      header=TRUE, stringsAsFactors = FALSE)
sharedStrainsSameSubjectByTimepointWide <-
  left_join(sharedStrainsSameSubjectByTimepointWide %>% mutate(subject=substr(sample,1,3)), 
            verifiedStrainTurnovers %>% dplyr::rename(subject=Subject, species_id=Species),
            by=c("subject", "species_id", "timeOfStrainTurnover"))
# Update the annotation of strain turnovers based on the SNP trajectory verification.
sharedStrainsSameSubjectByTimepointWide <- sharedStrainsSameSubjectByTimepointWide %>%
  mutate(strainTurnover=ifelse(!is.na(strainTurnoverVerifiedFilled), strainTurnoverVerifiedFilled, strainTurnover),
         timeOfStrainTurnover=ifelse(!is.na(strainTurnover) & strainTurnover, timeOfStrainTurnover, NA))
# Export the data that indicates at each timepoint whether or not a strain is shared with the initial timepoint
# and also whether it has experienced a strain turnover.
write.table(sharedStrainsSameSubjectByTimepointWide %>%
              dplyr::select(-strainTurnoverVerified, -strainTurnoverVerifiedFilled, -subject),
            paste0(OUTDIR, "sameSubject-sharedStrains-byTimepoint.txt"),
            row.names=FALSE, quote=FALSE)  

# Summarize the strain turnovers by species.
sharedStrainsSameSubjectTurnovers <- sharedStrainsSameSubjectByTimepointWide %>%
  mutate(subject=substr(sample,1,3), timepoint=as.integer(substr(sample,5,7))) %>%
  dplyr::select(subject, species_id, strainTurnover, timeOfStrainTurnover) %>%
  unique()

# For each species in each subject, determine if the initial strain
# is shared with various other timepoints (pre-abx, post-abx, follow-up, and final).
sharedStrainsSameSubjectBySpecies <- sharedStrainsSameSubject %>%
  filter(!is.na(consensusNA)) %>%
  group_by(species_id, subject1) %>%
  filter(sample1 %in% samplesInitial) %>%
  summarize(strainSharedInitialPreAbx=ifelse(sum(timepoint2<30)>0, 
                                      (sum(consensusNA[timepoint2<30])>0), NA),
            strainSharedInitialPostAbx=ifelse(sum(timepoint2>34 & timepoint2<75)>0, 
                                    (sum(consensusNA[timepoint2>34 & timepoint2<75])>0), NA),
            strainSharedInitialFollowup=ifelse(sum(timepoint2>75)>0, 
                                        (sum(consensusNA[timepoint2>75])>0), NA),
            strainSharedInitialPostAbxFollowup=ifelse(sum(timepoint2>34)>0, 
                                               (sum(consensusNA[timepoint2>34])>0), NA),
            strainSharedInitialFinal=ifelse(sum(sample2 %in% samplesFinal)>0, 
                                               (sum(consensusNA[sample2 %in% samplesFinal])>0), NA))
sharedStrainsSameSubjectBySpecies <- 
  left_join(sharedStrainsSameSubjectBySpecies %>% dplyr::rename(subject=subject1),
            sharedStrainsSameSubjectTurnovers,
            by=c("subject", "species_id"))
write.table(sharedStrainsSameSubjectBySpecies, 
            paste0(OUTDIR, "sameSubject-sharedStrains-subjectSummary.txt"),
            row.names=FALSE, quote=FALSE)

