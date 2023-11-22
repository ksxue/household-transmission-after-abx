# Import the annotated strain trajectory information.
# Summarize the total abundance of the disrupted and colonized species over time.
library(tidyverse)
library(RcppRoll)

# Import metadata about the study.
source("workflow/analysis/background/background.R")
# Set the output directory.
OUTDIR <- "workflow/analysis/annotateSpeciesTrajectoriesSharedStrains/out/"
LIMITOFDETECTION <- 1e-3

# Import the annotated species trajectories.
speciesAbundancesAnnotated <- 
  read.table(gzfile("workflow/analysis/annotateSpeciesTrajectoriesSharedStrains/out/speciesAbundances-trajectoriesSharingAnnotated-all.txt.gz"),
             header=TRUE, stringsAsFactors = FALSE)
# Import the summary of the annotated species trajectories.
speciesAbundancesAnnotatedSummary <-
  read.table("workflow/analysis/annotateSpeciesTrajectoriesSharedStrains/out/speciesTrajectoriesSharingSummary.txt",
             header=TRUE, stringsAsFactors = FALSE)

# Identify the new colonizers (including new strains and new species)..
colonizingStrainsSpecies <- speciesAbundancesAnnotatedSummary %>%
  filter(speciesStatusColonizationAllTimepoints!="speciesNotColonized" | strainTurnover) %>%
  mutate(timeOfColonization=ifelse(speciesStatusColonizationAllTimepoints!="speciesNotColonized",
                                   timeOfColonizationAllTimepoints, timeOfStrainTurnover),
         typeOfColonization=ifelse(speciesStatusColonizationAllTimepoints!="speciesNotColonized",
                                   "newSpecies", "newStrain")) %>%
  dplyr::select(subject, species_id, timeOfColonization, typeOfColonization)
# Export the list of new colonizers.
write.table(colonizingStrainsSpecies, paste0(OUTDIR, "colonizingStrainsSpecies.txt"),
            row.names=FALSE, quote = TRUE)

# Identify the resident species, i.e. species that were present before antibiotics.
# Calculate the appropriate background set of resident species.
# These have the same frequency properties as newly colonizing species
# but are not new colonizers.
# Import the species abundance data.
source("workflow/analysis/generateSpeciesAbundances/loadSpeciesAbundancesFiltered.R")
# Identify species that exceed a median relative abundance of 1e-3 before antibiotics
# but were detected at the beginning of the study.
# This contrasts with species that were not detected at the beginning of the study
# and later exceed a median relative abundance of 1e-3, which are dubbed new colonizers.
residentSpecies <- speciesAbundances %>%
  filter(timepoint<30) %>%
  group_by(subject, species_id) %>%
  mutate(window3Median=roll_median(relative_abundance, 3, align="left", fill=NA),
         window3Median=ifelse(is.na(window3Median), roll_median(relative_abundance, 2, align="left", fill=NA), window3Median),
         window3Median=ifelse(is.na(window3Median), relative_abundance, window3Median),
         timepointWindow3Middle=lag(timepoint, 1),
         timepointWindow3End=lag(timepoint, 2),
         window3Mean=roll_mean(relative_abundance, 3, align="left", fill=NA),
         initialWindow3Mean=window3Mean[sample %in% samplesInitial],
         maxWindow3MedianAll=max(window3Median)) %>%
  filter(initialWindow3Mean!=0 & maxWindow3MedianAll>1e-3) %>%
  group_by(subject, species_id) %>% 
  summarize(maxAbundance=max(relative_abundance),
            maxAbundanceMain=max(relative_abundance[timepoint<75]))
# Export the list of resident species.
write.table(residentSpecies, paste0(OUTDIR, "residentSpecies.txt"),
            row.names=FALSE, quote = TRUE)

# For each subject, summarize the total number of disrupted species
# based on calls from all timepoints and from the key timepoints.
numSpeciesDisruptedRecovered <- 
  rbind(speciesAbundancesAnnotatedSummary %>%
          group_by(subject, speciesStatusAbxAllTimepoints) %>%
          summarize(numSpecies=n()) %>% ungroup() %>% group_by(subject) %>%
          summarize(totalResident=sum(numSpecies[!is.na(speciesStatusAbxAllTimepoints)]),
                    disrupted=sum(numSpecies[speciesStatusAbxAllTimepoints %in% 
                                               c("speciesRecoveredFollowup", "speciesRecoveredMain", "speciesNotRecovered")]),
                    recoveredMain=sum(numSpecies[speciesStatusAbxAllTimepoints %in% c("speciesRecoveredMain")])) %>%
          mutate(timepoints="allTimepoints"),
        speciesAbundancesAnnotatedSummary %>%
          group_by(subject, speciesStatusAbxKTfollowup) %>%
          summarize(numSpecies=n()) %>% ungroup() %>% group_by(subject) %>%
          summarize(totalResident=sum(numSpecies[!is.na(speciesStatusAbxKTfollowup)]),
                    disrupted=sum(numSpecies[speciesStatusAbxKTfollowup %in% 
                                               c("speciesRecoveredFollowup", "speciesRecoveredMain", "speciesNotRecovered")]),
                    recoveredMain=sum(numSpecies[speciesStatusAbxKTfollowup %in% c("speciesRecoveredMain")])) %>%
          mutate(timepoints="keyTimepointsFollowup")) %>%
  pivot_longer(c("totalResident","disrupted","recoveredMain"), names_to="speciesTrajectoryAbx", values_to="numSpecies")
  
# For each subject, summarize the total number of colonized species
# based on calls from all timepoints and from the key timepoints.
numSpeciesColonized <- 
  rbind(speciesAbundancesAnnotatedSummary %>%
          group_by(subject, speciesStatusColonizationAllTimepoints) %>%
          summarize(numSpecies=n()) %>% ungroup() %>%
          complete(subject, speciesStatusColonizationAllTimepoints, fill=list(numSpecies=0)) %>%
          mutate(timepoints="allTimepoints") %>% 
          dplyr::rename(speciesStatusColonization=speciesStatusColonizationAllTimepoints),
        speciesAbundancesAnnotatedSummary %>%
          group_by(subject, speciesStatusColonizationKTfollowup) %>%
          summarize(numSpecies=n()) %>% ungroup() %>%
          complete(subject, speciesStatusColonizationKTfollowup, fill=list(numSpecies=0)) %>%
          mutate(timepoints="keyTimepointsFollowup") %>% 
          dplyr::rename(speciesStatusColonization=speciesStatusColonizationKTfollowup)) %>%
  filter(!is.na(speciesStatusColonization))

# For each subject, summarize the total number of strain turnovers.
# These strain turnovers are called based on all timepoints.
numStrainTurnovers <- speciesAbundancesAnnotatedSummary %>%
  filter(!is.na(strainTurnover)) %>%
  group_by(subject, strainTurnover) %>% 
  summarize(numSpecies=n()) %>% ungroup() %>%
  complete(subject, strainTurnover, fill=list(numSpecies=0)) %>%
  filter(strainTurnover) %>%
  mutate(timepoints="allTimepoints", strainTurnover="strainTurnover")

# Combine the species trajectory annotation summaries for disrupted, recovered, and colonized species.
speciesTrajectoriesSummary <-
  rbind(numSpeciesDisruptedRecovered %>%
          dplyr::rename(speciesTrajectory=speciesTrajectoryAbx),
        numSpeciesColonized %>%
          dplyr::rename(speciesTrajectory=speciesStatusColonization) %>%
          mutate(speciesTrajectory=
                   ifelse(speciesTrajectory=="speciesNotColonized", speciesTrajectory,
                   ifelse(speciesTrajectory=="pre-abx", "speciesColonizedPreAbx",
                   ifelse(speciesTrajectory=="abx", "speciesColonizedAbx",
                   ifelse(speciesTrajectory=="post-abx", "speciesColonizedPostAbx",
                   ifelse(speciesTrajectory=="followup", "speciesColonizedFollowup", NA)))))),
        numStrainTurnovers %>%
          dplyr::rename(speciesTrajectory=strainTurnover))
# Export species trajectory summary.
write.table(speciesTrajectoriesSummary, 
            paste0(OUTDIR, "speciesTrajectoriesPerSubject-summary.txt"),
            quote=FALSE, row.names=FALSE)

# For each subject, summarize the total relative abundance of the disrupted and colonized species over time.
# Also summarize the total relative abundance of each species with a new strain AFTER a strain turnover occurs
# (that is, do not count the initial strain in this relative abundance total).
# Use the annotations from all sequenced timepoints.
totalAbundanceDisruptedColonized <- speciesAbundancesAnnotated %>%
  group_by(hh, subject, sample, timepoint) %>%
  summarize(disrupted=
              sum(relative_abundance[speciesStatusAbxAllTimepoints %in%
                                    c("speciesRecoveredFollowup", "speciesRecoveredMain", "speciesNotRecovered")]),
            colonized=
             sum(relative_abundance[!is.na(speciesStatusColonizationAllTimepoints) & speciesStatusColonizationAllTimepoints!="speciesNotColonized"]),
            strainTurnover=
              sum(relative_abundance[!is.na(strainTurnover) & strainTurnover & timepoint>=timeOfStrainTurnover]),
            disruptedNotRecovered=
              sum(relative_abundance[speciesStatusAbxAllTimepoints %in% c("speciesNotRecovered", "speciesRecoveredFollowup")])) %>%
  pivot_longer(c("disrupted","colonized","strainTurnover","disruptedNotRecovered"), 
               names_to="speciesTrajectory", values_to="totalAbundance")
# Export the summary of the total abundances of the disrupted and colonizing species,
# as well as the species that undergo strain turnovers.
write.table(totalAbundanceDisruptedColonized,
            paste0(OUTDIR, "speciesTrajectoriesPerSubject-totalAbundance-allTimepoints.txt"),
            quote=FALSE, row.names=FALSE)

# For each subject, summarize the total relative abundance of the disrupted and colonized species over time.
# Also summarize the total relative abundance of each species with a new strain AFTER a strain turnover occurs
# (that is, do not count the initial strain in this relative abundance total).
# Use the annotations using only the main timepoints.
totalAbundanceDisruptedColonizedKTfollowup <- speciesAbundancesAnnotated %>%
  group_by(hh, subject, sample, timepoint) %>%
  summarize(disrupted=
              sum(relative_abundance[speciesStatusAbxKTfollowup %in%
                                       c("speciesRecoveredFollowup", "speciesRecoveredMain", "speciesNotRecovered")]),
            colonized=sum(relative_abundance[!is.na(speciesStatusColonizationKTfollowup) &
                                               speciesStatusColonizationKTfollowup!="speciesNotColonized"]),
            strainTurnover=
              sum(relative_abundance[!is.na(strainTurnover) & strainTurnover & timepoint>=timeOfStrainTurnover]),
            disruptedNotRecovered=
              sum(relative_abundance[speciesStatusAbxAllTimepoints %in% c("speciesNotRecovered", "speciesRecoveredFollowup")])) %>%
  pivot_longer(c("disrupted","colonized","strainTurnover","disruptedNotRecovered"), 
               names_to="speciesTrajectory", values_to="totalAbundance")
# Export the summary of the total abundances of the disrupted and colonizing species,
# as well as the species that undergo strain turnovers.
write.table(totalAbundanceDisruptedColonizedKTfollowup,
            paste0(OUTDIR, "speciesTrajectoriesPerSubject-totalAbundance-keyTimepointsFollowup.txt"),
            quote=FALSE, row.names=FALSE)
