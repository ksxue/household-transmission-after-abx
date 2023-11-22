# This script runs all analyses comparing the species abundances
# of samples at the initial timepoint of the study.
# It requires the output of generateSpeciesAbundances to be present.

# This script generates a list of all pairs of samples collected at the initial timepoint,
# and it calculates the number of species and the relative abundance of species
# shared above a limit of detection for each pair of samples.
source("workflow/analysis/compareSpeciesAbundances/compareSpeciesAbundances.R")
rm(list=ls())

# Touch a file to verify that the script finished running.
file.create("workflow/analysis/compareSpeciesAbundances/done.txt")
