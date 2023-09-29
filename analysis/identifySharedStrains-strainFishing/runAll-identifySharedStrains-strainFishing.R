# This script runs all analyses associated with identifying shared strains
# using strain fishing data.
# The script requires the output of the snakemake performStrainFishing rule
# to be present.

# Iterate through all analyzed species and aggregate the strain fishing calls.
# Output this set of strain fishing calls.
source("workflow/analysis/identifySharedStrains-strainFishing/generateStrainFishingDataframe.R")
rm(list=ls())

# Identify shared strains based on the strain fishing calls.
# Output the shared strains, strain cousins, and number of private SNPs per bait.
source("workflow/analysis/identifySharedStrains-strainFishing/identifySharedStrains-strainFishing.R")
rm(list=ls())

# Summarize the strain sharing calls based on strain fishing.
# Combine the strain sharing calls from both "directions" where applicable.
# Output the strain sharing calls at the initial timepoint and later key timepoints.
source("workflow/analysis/identifySharedStrains-strainFishing/summarizeSharedStrains-strainFishing.R")
rm(list=ls())

# Touch a file to verify that the script finished running.
file.create("workflow/analysis/identifySharedStrains-strainFishing/done.txt")
