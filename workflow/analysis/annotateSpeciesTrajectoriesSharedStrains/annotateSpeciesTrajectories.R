library(tidyverse)
library(RcppRoll)

# Import metadata about the study.
source("workflow/analysis/background/background.R")
# Set the output directory.
OUTDIR <- "workflow/analysis/annotateSpeciesTrajectoriesSharedStrains/out/"
LIMITOFDETECTION <- 1e-3


# Import species abundances and annotate species trajectories -------------

# Import the species abundance data.
source("workflow/analysis/generateSpeciesAbundances/loadSpeciesAbundancesFiltered.R")

# Write a function that, given a dataframe of species abundances,
# calculates various mean and median abundances for different parts of the study.
annotateSpeciesAbundances <- function(ispeciesAbundances){
  # Calculate the pre- and post-abx median and minimum abundances.
  speciesAbundancesAnnotated <- ispeciesAbundances %>%
    group_by(subject, species_id) %>%
    mutate(medianPreAbxAbundance=median(relative_abundance[timepoint<30]),
           medianPostAbxAbundance=median(relative_abundance[timepoint>34 & timepoint<75]),
           minPostAbxAbundance=min(relative_abundance[timepoint>29 & timepoint<75]),
           timeOfMinPostAbxAbundance=max(timepoint[relative_abundance==minPostAbxAbundance & timepoint<75]))
  # Set the median and min abundances to 1e-6 if they are 0 to prevent problems with ratios.
  speciesAbundancesAnnotated <- speciesAbundancesAnnotated %>%
    mutate(medianPreAbxAbundance=ifelse(medianPreAbxAbundance==0,1e-6,medianPreAbxAbundance),
           medianPostAbxAbundance=ifelse(medianPostAbxAbundance==0,1e-6,medianPostAbxAbundance),
           minPostAbxAbundance=ifelse(minPostAbxAbundance==0,1e-6,minPostAbxAbundance))
  
  # Calculate a rolling median of the species abundances in a 3-sample window.
  # Implement a partial window median at the end of the sampling timeframe
  speciesAbundancesAnnotated <- speciesAbundancesAnnotated %>%
    group_by(subject, species_id) %>%
    mutate(window3Median=roll_median(relative_abundance, 3, align="left", fill=NA),
           window3Median=ifelse(is.na(window3Median), roll_median(relative_abundance, 2, align="left", fill=NA), window3Median),
           window3Median=ifelse(is.na(window3Median), relative_abundance, window3Median),
           window3Mean=roll_mean(relative_abundance, 3, align="left", fill=NA),
           timepointWindow3Middle=lag(timepoint, 1),
           timepointWindow3End=lag(timepoint, 2)) %>%
    mutate(maxPreAbxWindow3Median=max(window3Median[timepoint<30]),
           minPostAbxWindow3Median=min(window3Median[timepoint>34 & timepoint<75]),
           maxPostAbxWindow3Median=max(window3Median[timepoint>34 & timepoint<75]),
           maxPostAbxAllWindow3Median=max(window3Median[timepoint>34]),
           initialWindow3Median=window3Median[sample %in% samplesInitial],
           initialWindow3Mean=window3Mean[sample %in% samplesInitial],
           maxWindow3MedianMain=max(window3Median[timepoint<75]),
           maxWindow3MedianAll=max(window3Median))
  
  # Annotate species based on disruption, recovery, time to recovery,
  # colonization, and time of colonization.
  speciesAbundancesAnnotated <- speciesAbundancesAnnotated %>%
    group_by(subject, species_id) %>%
    # Require that species have a median pre-abx abundance above the limit of detection
    # in order to reliably detect whether the species are disrupted.
    # It is difficult to detect whether low-abundance species have truly been disrupted.
    mutate(speciesDisrupted=ifelse(maxPreAbxWindow3Median>LIMITOFDETECTION, minPostAbxAbundance/medianPreAbxAbundance<10^-1.5, NA),
           speciesRecovered=(speciesDisrupted & 
                               max(maxPostAbxAllWindow3Median/medianPreAbxAbundance>10^-1)),
           # When calculating time of recovery or colonization, take the later of the
           # time window or the timepoint at which relative abundance exceeds the limit of detection
           # to account for "edge effects" in which the window spans a wide range of time
           # and the left alignment of the window may give an incorrect idea of when an event occurred.
           timeOfRecovery=ifelse(!speciesRecovered, NA,
                                 max(min(timepoint[timepoint>34 & window3Median/medianPreAbxAbundance>10^-1]),
                                     min(timepoint[timepoint>34 & relative_abundance/medianPreAbxAbundance>10^-1]))),
           speciesColonized=(initialWindow3Mean==0 & maxWindow3MedianAll>LIMITOFDETECTION),
           timeOfColonization=ifelse(!speciesColonized, NA,
                                     max(min(timepoint[window3Median>LIMITOFDETECTION]),
                                         min(timepoint[relative_abundance>LIMITOFDETECTION]))))
  # Annotate whether a species colonizes transiently
  # or becomes a persistent member of the gut microbiome.
  # Determine this based on whether the species is detected at the final timepoint.
  speciesAbundancesAnnotated <- speciesAbundancesAnnotated %>%
    group_by(subject, species_id) %>%
    mutate(speciesColonizedTransientlyMain=
             ifelse(!(speciesColonized & timeOfColonization<75), NA,
                    ifelse(relative_abundance[sample %in% samplesFinal]<LIMITOFDETECTION, TRUE, FALSE)),
           speciesColonizedTransientlyAll=
             ifelse(!(speciesColonized), NA,
                    ifelse(relative_abundance[sample %in% samplesLastSequencedX]<LIMITOFDETECTION, TRUE, FALSE)))
  
  # Retain only species whose median over a 3-sample window exceeds the limit of detection.
  speciesAbundancesAnnotated <- speciesAbundancesAnnotated %>%
    filter(maxWindow3MedianAll>LIMITOFDETECTION)
  
  # Annotate species behaviors as disrupted, recovered, and colonized,
  # and also summarize the timing of each of these behaviors.
  speciesAbundancesAnnotated <- speciesAbundancesAnnotated %>%
    mutate(speciesStatusAbx=ifelse(!speciesDisrupted, "speciesNotDisrupted",
                                   ifelse(!speciesRecovered, "speciesNotRecovered",
                                          ifelse(timeOfRecovery<75, "speciesRecoveredMain", "speciesRecoveredFollowup")))) %>%
    mutate(speciesStatusColonization=ifelse(!speciesColonized, "speciesNotColonized",
                                            ifelse(timeOfColonization<30, "pre-abx",
                                                   ifelse(timeOfColonization<36, "abx",
                                                          ifelse(timeOfColonization<75, "post-abx", "followup"))))) %>%
    mutate(speciesStatus=ifelse(is.na(speciesStatusAbx) & speciesStatusColonization=="speciesNotColonized", "none",
                                ifelse(is.na(speciesStatusAbx), paste0("speciesColonized-", speciesStatusColonization), speciesStatusAbx)))
  
  return(speciesAbundancesAnnotated)
}

# Retain all samples for subjects in study arm X.
# Focus only on species that exceed the limit of detection.
speciesAbundancesXall <- speciesAbundances %>%
  group_by(subject, species_id) %>%
  filter(max(relative_abundance)>LIMITOFDETECTION, subject %in% subjectsX)
speciesAbundancesAnnotatedXall <- annotateSpeciesAbundances(speciesAbundancesXall)
# Remove columns used to annotate species trajectories.
speciesAbundancesAnnotatedXall <- speciesAbundancesAnnotatedXall %>%
  dplyr::select(-medianPreAbxAbundance, -medianPostAbxAbundance, -minPostAbxAbundance,
                -timeOfMinPostAbxAbundance, -window3Median, -window3Mean, -maxPreAbxWindow3Median,
                -minPostAbxWindow3Median, -maxPostAbxWindow3Median, -maxPostAbxAllWindow3Median,
                -initialWindow3Median, -initialWindow3Mean, -maxWindow3MedianMain, -maxWindow3MedianAll)
# Summarize the overall behavior of each species.
speciesAbundancesAnnotatedSummaryXall <- speciesAbundancesAnnotatedXall %>%
  dplyr::select(subject, species_id,
                kingdom, phylum, class, order, family, genus, species,
                speciesDisrupted, speciesRecovered, speciesColonized,
                speciesStatusAbx, speciesStatusColonization, 
                timeOfRecovery, timeOfColonization, 
                speciesColonizedTransientlyMain, speciesColonizedTransientlyAll) %>%
  unique()
# Export the full abundance info annotated with species trajectories.
write.table(speciesAbundancesAnnotatedXall, 
            gzfile(paste0(OUTDIR, "speciesAbundances-trajectoriesAnnotated-all.txt.gz")),
            row.names=FALSE)
# Export the summary of the species trajectory behavior.
write.table(speciesAbundancesAnnotatedSummaryXall, 
            paste0(OUTDIR, "speciesTrajectoriesSummary-all.txt"),
            row.names=FALSE)

# Retain only key timepoints and follow-up samples for subjects in study arm X.
# Focus only on species that exceed the limit of detection.
speciesAbundancesXKTfollowup <- speciesAbundances %>%
  group_by(subject, species_id) %>%
  filter(max(relative_abundance)>LIMITOFDETECTION, subject %in% subjectsX,
         sample %in% c(samplesKeyTimepoints, samplesXfollowup))
speciesAbundancesAnnotatedXKTfollowup <- annotateSpeciesAbundances(speciesAbundancesXKTfollowup)
# Summarize the overall behavior of each species.
speciesAbundancesAnnotatedSummaryXKTfollowup <- speciesAbundancesAnnotatedXKTfollowup %>%
  dplyr::select(subject, species_id,
                kingdom, phylum, class, order, family, genus, species,
                speciesDisrupted, speciesRecovered, speciesColonized,
                speciesStatusAbx, speciesStatusColonization, 
                timeOfRecovery, timeOfColonization) %>%
  unique()
# Export the full abundance info annotated with species trajectories.
write.table(speciesAbundancesAnnotatedXKTfollowup, 
            gzfile(paste0(OUTDIR, "speciesAbundances-trajectoriesAnnotated-keyTimepointsFollowup.txt.gz")),
            row.names=FALSE)
# Export the summary of the species trajectory behavior.
write.table(speciesAbundancesAnnotatedSummaryXKTfollowup, 
            paste0(OUTDIR, "speciesTrajectoriesSummary-keyTimepointsFollowup.txt"),
            row.names=FALSE)
