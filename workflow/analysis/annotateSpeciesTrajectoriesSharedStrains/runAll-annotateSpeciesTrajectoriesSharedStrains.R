# This script annotates the species relative abundance trajectories.
source("workflow/analysis/annotateSpeciesTrajectoriesSharedStrains/annotateSpeciesTrajectories.R")
rm(list=ls())

# This script annotates the species sharing between subjects.
source("workflow/analysis/annotateSpeciesTrajectoriesSharedStrains/annotateSpeciesSharing.R")
rm(list=ls())

# This script integrates calls of species trajectory behaviors (disrupted, recovered, colonized)
# with calls of species sharing and calls of strain sharing and strain turnovers.
# It also makes use of additional SNP data to verify candidate strain turnovers.
source("workflow/analysis/annotateSpeciesTrajectoriesSharedStrains/integrateSpeciesTrajectoriesStrainSharing.R")
rm(list=ls())

# This script takes in the calls of species trajectories and strain sharing
# and summarizes the total proportion of the microbiome comprised of
# disrupted, recovered, and colonizing species and strains.
source("workflow/analysis/annotateSpeciesTrajectoriesSharedStrains/summarizeSpeciesTrajectoriesStrainSharing.R")
rm(list=ls())

# Touch a file to verify that the script finished running.
file.create("workflow/analysis/annotateSpeciesTrajectoriesSharedStrains/done.txt")
