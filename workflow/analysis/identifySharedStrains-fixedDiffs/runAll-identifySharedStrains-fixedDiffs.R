# This script runs all analyses associated with identifying shared strains
# based on the number of fixed differences between populations.
# This script requires the output of identifyContaminatedSamples
# to be present.

# Load data about the number of fixed differences between populations.
# Filter out samples with evidence of contamination or that are on the sample blacklist.
# Set empirical filters based on the between-household distance distribution for each species
# and annotate each pair of populations as shared or not.
source("workflow/analysis/identifySharedStrains-fixedDiffs/identifySharedStrains-fixedDiffs.R")
rm(list=ls())

# Load the annotated dataset of fixed differences between pairs of samples within households.
# Summarize and output the number of shared strains at the initial timepoint and other key timepoints.
source("workflow/analysis/identifySharedStrains-fixedDiffs/summarizeSharedStrains-fixedDiffs.R")
rm(list=ls())

# Touch a file to verify that the script finished running.
file.create("workflow/analysis/identifySharedStrains-fixedDiffs/done.txt")