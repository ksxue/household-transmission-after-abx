library(tidyverse)

# Import background metadata and plot themes.
source("workflow/analysis/background/background.R")
source("workflow/analysis/plotDefaults.R")
LIMITOFDETECTION <- 1e-5
OUTDIR <- "workflow/analysis/fitColonizationTrajectories/out/"


# Parse the species and strains to identify colonization and recov --------

# Import annotated species trajectories.
dataSpeciesTrajectories <- 
  read.table("workflow/analysis/annotateSpeciesTrajectoriesSharedStrains/out/speciesAbundances-trajectoriesSharingAnnotated-all.txt.gz",
             header=TRUE, stringsAsFactors = FALSE)

# Import annotated species trajectory summaries.
dataSpeciesTrajectoriesSummary <-
  read.table("workflow/analysis/annotateSpeciesTrajectoriesSharedStrains/out/speciesTrajectoriesSharingSummary.txt",
             header=TRUE, stringsAsFactors = FALSE)

# Identify the species that have recovered, colonized, or experienced a species turnover.
dataSpeciesSummaryRecoveryColonization <- dataSpeciesTrajectoriesSummary %>%
  filter(!is.na(timeOfRecoveryAllTimepoints) | !is.na(timeOfColonizationAllTimepoints) | !is.na(timeOfStrainTurnover)) %>%
  dplyr::select(subject, species_id, speciesStatusAbxAllTimepoints, speciesStatusColonizationAllTimepoints,
                timeOfRecoveryAllTimepoints, timeOfColonizationAllTimepoints, 
                strainTurnover, timeOfStrainTurnover, sharingHhPreAbx, sharingHhPostAbx, sharingHhFollowup,
                sharingSubjectPreAbx, sharingSubjectPostAbx, sharingSubjectFollowup)
# Infer the time at which this event occurred.
dataSpeciesSummaryRecoveryColonization <- dataSpeciesSummaryRecoveryColonization %>%
  mutate(timeDetected=case_when(
    !is.na(timeOfRecoveryAllTimepoints) ~ timeOfRecoveryAllTimepoints,
    !is.na(timeOfColonizationAllTimepoints) ~ timeOfColonizationAllTimepoints,
    !is.na(timeOfStrainTurnover) ~ timeOfStrainTurnover
  ))
# Annotate the strain sharing status for each event.
dataSpeciesSummaryRecoveryColonization <- dataSpeciesSummaryRecoveryColonization %>%
  mutate(strainType=case_when(
    !is.na(strainTurnover) & strainTurnover ~ "new strain",
    !is.na(timeOfColonizationAllTimepoints) ~ "new species",
    timeDetected<75 & sharingSubjectPostAbx=="strainsShared" ~ "same strain",
    timeDetected>75 & sharingSubjectFollowup=="strainsShared" ~ "same strain",
    timeDetected<75 & sharingSubjectPostAbx=="strainsUnknown" ~ "unknown strain",
    timeDetected>75 & sharingSubjectFollowup=="strainsUnknown" ~ "unknown strain"
  ))

# Extract the full species trajectories for these recovery and colonization events.
dataSpeciesTrajectoriesRecoveryColonization <- 
  inner_join(dataSpeciesTrajectories,
             dataSpeciesSummaryRecoveryColonization %>% 
               dplyr::select(subject, species_id, timeDetected, strainType) %>%
               dplyr::rename(timeOfRecoveryColonizationTurnover=timeDetected,
                             typeOfRecoveryColonizationTurnover=strainType))
# Simplify the columns attached to the full species trajectories.
dataSpeciesTrajectoriesRecoveryColonization <- dataSpeciesTrajectoriesRecoveryColonization %>%
  dplyr::select(sample, subject, hh, timepoint, species_id, relative_abundance, totalReads,
               speciesColonizedTransientlyMain, speciesColonizedTransientlyAll,
               timeOfRecoveryColonizationTurnover, typeOfRecoveryColonizationTurnover)
# Remove unnecessary columns from these trajectories and export the species relative abundances.
write.table(dataSpeciesTrajectoriesRecoveryColonization,
            paste0(OUTDIR, "speciesRecoveryColonizationStrainTurnovers-trajectories.txt"),
            row.names=FALSE, quote=FALSE, sep="\t")


# Import the trajectory fitting parameters --------------------------------

# Import trajectory fitting parameters.
dataTrajectories <- read.table("workflow/analysis/fitColonizationTrajectories/out/trajectoryFits/output_summary/trajectory_fits.txt",
                               header=TRUE, stringsAsFactors = FALSE)
# Filter out species that don't have trajectory fits.
dataTrajectoriesFit <- dataTrajectories %>%
  filter(Reason=="SUCCESS")
dataTrajectoriesFit <- dataTrajectoriesFit %>%
  mutate(tstarTimepoint=
           ifelse(abs(tstar-closest_t_before_tstar)<abs(tstar-closest_t_after_star),
                  closest_t_before_tstar, closest_t_after_star))

# Append the trajectory fit parameters to the full species abundance trajectories.
dataSpeciesAbundancesTrajectoryFit <-
  left_join(dataSpeciesTrajectoriesRecoveryColonization,
            dataTrajectoriesFit %>% dplyr::rename(subject=Subject, species_id=Species),
            by=c("subject","species_id"))

# For each species with a trajectory fit,
# calculate the number of samples that you have to go back from the tstarTimepoint
# (inferred time of saturation) for the species to go below the limit of detection.
# In other words, calculate the number of timepoints required to go from below the limit
# of detection to the tstarTimepoint.
dataSpeciesAbundancesTrajectoryFit <- dataSpeciesAbundancesTrajectoryFit %>%
  group_by(subject, species_id) %>%
  arrange(subject, species_id, timepoint) %>%
  mutate(timepointIndex=row_number())
# Calculate the last sample before tstar at which the sample is below the limit of detection.
dataSpeciesAbundancesTrajectoryFit <- dataSpeciesAbundancesTrajectoryFit %>%
  filter(Reason=="SUCCESS") %>%
  group_by(subject, species_id) %>%
  mutate(timepointIndexTstar=timepointIndex[timepoint==tstarTimepoint],
         lastUndetectableTimepointIndexBeforeTstar=
           max(timepointIndex[timepoint<tstarTimepoint & relative_abundance<LIMITOFDETECTION]),
         lastUndetectableTimepointBeforeTstar=
           max(timepoint[timepoint<tstarTimepoint & relative_abundance<LIMITOFDETECTION]),
         lastTimepointBeforeTstar=
           timepoint[timepointIndex==timepointIndexTstar-1])
# Calculate the number of timepoints it takes to increase from undetectable levels to K.
dataSpeciesAbundancesTrajectoryFit <- dataSpeciesAbundancesTrajectoryFit %>%
  mutate(numTimepointsToReachK=timepointIndexTstar-lastUndetectableTimepointIndexBeforeTstar,
         timeIntervalToReachK=tstarTimepoint-lastUndetectableTimepointBeforeTstar)
# Summarize the trajectory fits and waiting time for each species.
dataSpeciesAbundancesTrajectoryFitSummary <- dataSpeciesAbundancesTrajectoryFit %>%
  ungroup() %>%
  dplyr::select(subject, species_id, timeOfRecoveryColonizationTurnover, typeOfRecoveryColonizationTurnover,
                speciesColonizedTransientlyMain, speciesColonizedTransientlyAll,
                r, K, tstar, closest_t_before_tstar, closest_t_after_star,
                inferred_tau, inferred_f0, Reason, 
                tstarTimepoint, timepointIndexTstar, lastTimepointBeforeTstar,
                lastUndetectableTimepointIndexBeforeTstar, lastUndetectableTimepointBeforeTstar,
                numTimepointsToReachK, timeIntervalToReachK) %>%
  unique()

# Export the annotated species trajectories.
write.table(dataSpeciesAbundancesTrajectoryFit, paste0(OUTDIR, "speciesTrajectoriesFit.txt"),
            row.names=FALSE, quote=FALSE, sep="\t")
write.table(dataSpeciesAbundancesTrajectoryFitSummary, 
            paste0(OUTDIR, "speciesTrajectoriesFit-summary.txt"),
            row.names=FALSE, quote=FALSE, sep="\t")
