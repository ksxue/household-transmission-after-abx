# This script runs all analyses associated with identifying contaminated samples.
# The script requires the output of the snakemake calculateFixedDiffs rule
# to be present.

# Aggregate all fixed differences output files.
source("workflow/analysis/identifyContaminatedSamples/generateFixedDiffsDataframe.R")
rm(list=ls())

# Identify contaminated samples.
# This script also calls the loadFixedDiffs.R file.
source("workflow/analysis/identifyContaminatedSamples/identifyContaminatedSamples.R")
rm(list=ls())

# Touch a file to verify that the script finished running.
file.create("workflow/analysis/identifyContaminatedSamples/done.txt")