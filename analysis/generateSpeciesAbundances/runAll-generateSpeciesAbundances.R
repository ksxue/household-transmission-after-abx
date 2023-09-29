# This script imports calculations of species abundances from MIDAS,
# combines the counts of single-copy genes with the taxonomy information,
# and outputs a dataframe summarizing species abundances.
# The script requires the background module to have completed running.

# Import species abundance calls from MIDAS. Merge with species taxonomy
# and output a consolidated set of species abundance calls.
source("workflow/analysis/generateSpeciesAbundances/generateSpeciesAbundances.R")
rm(list=ls())

# Generate a palette for coloring the most common bacterial families.
source("workflow/analysis/generateSpeciesAbundances/generateSpeciesPalettes.R")
rm(list=ls())

# Touch a file to verify that the script finished running.
file.create("workflow/analysis/compareSharedStrainCalls/done.txt")
