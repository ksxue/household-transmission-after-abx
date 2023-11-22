# This script takes in the inferred timing of colonization events
# based on the trajectory fitting data.
# It analyzes the clustering of colonization events during follow-up sampling
# and simulates colonization events under a null model of uniformly distributed
# colonization events to test whether the maximum clustering of colonization events
# in one subject is more than would be expected due to chance.
library(tidyverse)

# Import the study metadata.
source("workflow/analysis/background/background.R")
# Set the output directory.
OUTDIR <- "workflow/analysis/analyzeClusteringOfColonizationEvents/out/"
# Set the number of simulated sets of colonization events.
NUMSIMULATIONS <- 1000


# Import the classification of subject responses.
subjectResponses <- read.table("workflow/analysis/classifySubjectResponses/out/subjectResponses.txt",
                               header=TRUE, stringsAsFactors = FALSE, sep="\t")

# Import species trajectory fitting summaries.
dataSpeciesTrajectoryFitSummary <-
  read.table("workflow/analysis/fitColonizationTrajectories/out/speciesTrajectoriesFit-summary.txt",
             header=TRUE, stringsAsFactors = TRUE, sep="\t")
# Extract the list of colonizers and strain turnovers.
# Exclude unknown strains, which do not have enough depth of coverage to determine if they are new or old strains.
dataColonizationTurnoverTiming <- dataSpeciesTrajectoryFitSummary %>%
  filter(typeOfRecoveryColonizationTurnover %in% c("new strain", "new species") &
           tstarTimepoint>75) %>%
  dplyr::select(subject, species_id, typeOfRecoveryColonizationTurnover, tstarTimepoint) %>%
  dplyr::rename(typeOfEvent=typeOfRecoveryColonizationTurnover,
                timingOfEvent=tstarTimepoint)

# Write a function that takes in the name of a subject,
# extracts the list of sampled timepoints and the list of colonization/turnover event timing,
# and uses permutation tests to determine the maximum number of events at one timepoint.
# Use a null model in which colonization events are distributed uniformly
# during follow-up sampling.
randomizeEventTiming <- function(isubject){
  # Extract the list of follow-up sampling timepoints.
  samplesFollowup <- samplesRaw %>%
    filter(subject==isubject, timepoint>75) %>% arrange(timepoint) %>% pull(timepoint)
  # Extract the list of colonization and strain turnover events during follow-up sampling.
  colonizationEventsFollowup <- dataColonizationTurnoverTiming %>%
    filter(subject==isubject, timingOfEvent>75) %>% arrange(timingOfEvent)
  
  # Calculate the main parameters for running the permutation tests.
  # This includes the duration of follow-up sampling, the number of colonization events,
  # and the maximum number of colonization events detected at one timepoint.
  durationOfFollowupSampling <- 
    max(samplesRaw %>% filter(subject==isubject) %>% pull(timepoint)) -
    max(samplesRaw %>% filter(subject==isubject, sample %in% samplesXmain) %>% pull(timepoint))
  numColonizationEvents <- nrow(colonizationEventsFollowup)
  maxNumColonizationEventsOneTimepoint <- unique(dataColonizationTurnoverTiming %>%
                                                   filter(subject==isubject, timingOfEvent>75) %>%
                                                   group_by(timingOfEvent) %>% summarize(numEvents=n()) %>% ungroup() %>%
                                                   filter(numEvents==max(numEvents)) %>% pull(numEvents))
  # Identify the timepoint with a cohort of colonization events of this size
  # that has the shortest time interval preceding it.
  timeOfMaxColonizationEvents <- dataColonizationTurnoverTiming %>%
    filter(subject==isubject, timingOfEvent>75) %>%
    group_by(timingOfEvent) %>% summarize(numEvents=n()) %>% ungroup() %>%
    arrange(timingOfEvent) %>%
    mutate(timeToPreviousSample=timingOfEvent-lag(timingOfEvent, default=0)) %>%
    filter(numEvents==max(numEvents)) %>%
    filter(timeToPreviousSample==min(timeToPreviousSample)) %>%
    pull(timingOfEvent)
  
  # Set seed to generate consistent simulations.
  set.seed(0)
  # Simulate 1000 distributions of colonization events and determine the number
  # of colonization events detected at the timepoint that has the highest number of
  # real colonization events.
  simulatedNumbersOfColonizationEvents <- foreach(i=seq(1:NUMSIMULATIONS), .combine="c") %do% {
    # Draw a set of uniformly distributed random variables.
    # Convert the variables to colonization timing during follow-up sampling.
    simulatedColonizationEvents <- sort(runif(numColonizationEvents)*durationOfFollowupSampling+
                                          max(samplesRaw %>% filter(subject==isubject, sample %in% samplesXmain) %>% pull(timepoint)))
    
    # Calculate the number of events that would be detected at each follow-up sampling
    # timepoint given the timing of follow-up sampling points.
    totalSimulatedEventsPerFollowupTimepoint <- 
      sapply(sort(samplesFollowup), function(x) length(which(simulatedColonizationEvents<x)))
    # Given the list of cumulative events, also identify the number of new events
    # detected at each follow-up sampling timepoint.
    newSimulatedEventsPerFollowupTimepoint <- 
      totalSimulatedEventsPerFollowupTimepoint - c(0,head(totalSimulatedEventsPerFollowupTimepoint,-1))
    names(newSimulatedEventsPerFollowupTimepoint) <- as.character(samplesFollowup)
    # Identify the list of simulated events at the timepoint that has the
    # highest number of colonization events in the real follow-up sampling data.
    numSimulatedEventsDetectedAtTimepointOfMaxRealColonizationEvents <- 
      newSimulatedEventsPerFollowupTimepoint[as.character(timeOfMaxColonizationEvents)]
  }
  
  # Convert vector of the number of simulated events to a dataframe.
  return(as.data.frame(simulatedNumbersOfColonizationEvents) %>%
           mutate(subject=isubject, timeOfColonizationEvents=timeOfMaxColonizationEvents,
                  actualNumberOfColonizationEvents=maxNumColonizationEventsOneTimepoint))
  
}

dataSimulationsAllSubjects <- 
  foreach(isubject=(subjectResponses %>% filter(subjectResponse=="lasting response") %>% pull(subject)), 
          .combine="rbind") %do% {
            randomizeEventTiming(isubject)
          }

# Export simulations and actual colonization timing for subjects with lasting responses.
write.table(dataSimulationsAllSubjects, paste0(OUTDIR, "simulations-colonizationCohorts.txt"),
            quote=FALSE, row.names=FALSE)

