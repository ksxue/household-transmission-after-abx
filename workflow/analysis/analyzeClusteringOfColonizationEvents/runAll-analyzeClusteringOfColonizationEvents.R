# This module takes the list of observed colonization events
# and performs permutation tests to determine if there are more events observed at a timepoint
# than predicted under uniform expectations.
# Perform permutation tests.
source("workflow/analysis/analyzeClusteringOfColonizationEvents/analyzeClusteringofColonizationEvents.R")
rm(list=ls())
