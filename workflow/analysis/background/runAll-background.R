# This script runs all component scripts in the "background" folder.
# Note that generateSampleBlacklist.R makes use of background.R.

source("workflow/analysis/background/generateSampleBlacklist.R")
rm(list=ls())

# Touch a file to verify that the script finished running.
file.create("workflow/analysis/background/done.txt")
