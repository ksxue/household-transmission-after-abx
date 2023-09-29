# This module takes the HMP species abundances as processed by Good and Rosenberg,
# calculates the prevalence and abundance of each gut species,
# and calculates the HMP prevalence and abundance of resident species and new colonizers.

# Calculate the prevalence and abundance of species in the HMP.
source("workflow/analysis/calculateHMPspeciesPrevalenceAbundance/calculateHMPspeciesPrevalenceAbundance.R")
rm(list=ls())


# Compare the HMP prevalence and abundance of resident species and new colonizers. 
source("workflow/analysis/calculateHMPspeciesPrevalenceAbundance/compareResidentSpeciesNewColonizers.R")
rm(list=ls())