# This module brings in the species abundances calculated by MIDAS
# and annotates the relative abundances in sample pairs based on strain sharing.

# This script calculates the number and proportion of species and strains shared
# between pairs of samples from matched key timepoints from cohabiting partners.
source("workflow/analysis/integrateSpeciesAbundancesSharedStrains/integrateSpeciesAbundancesSharedStrains.R")
rm(list=ls())

# This script calculates the number and proportion of species and strains shared
# between pairs of samples from the initial timepoint from cohabiting and non-cohabiting partners.
source("workflow/analysis/integrateSpeciesAbundancesSharedStrains/integrateSpeciesAbundancesSharedStrains-allInitialSamples.R")
rm(list=ls())

# This script calculates the number and proportion of species and strain shared
# between pairs of samples from the same subject over time.
source("workflow/analysis/integrateSpeciesAbundancesSharedStrains/integrateSpeciesAbundancesSharedStrains-initialFinalSameSubject.R")
rm(list=ls())
