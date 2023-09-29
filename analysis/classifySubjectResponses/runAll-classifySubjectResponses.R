# This script uses the JSD from the initial timepoint to calculate maximum and final JSD
# and classify antibiotic-taking subjects based on antibiotic response.
source("workflow/analysis/classifySubjectResponses/classifySubjectResponses.R")
rm(list=ls())

# Touch a file to verify that the script finished running.
file.create("workflow/analysis/classifySubjectResponses/done.txt")
