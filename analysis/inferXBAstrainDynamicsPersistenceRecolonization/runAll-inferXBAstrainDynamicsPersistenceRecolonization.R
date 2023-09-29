# Perform SNP analyses to determine whether the increase in abundance
# of XBA strains after antibiotics is due to strain persistence or recolonization.

# Analyze the SNP cohorts in XBA before and after antibiotics.
# Identify SNPs that distinguish the XBA and XBB strains,
# SNPs that fix after antibiotics, and also SNPs that distinguish
# co-colonizing strains.
source("workflow/analysis/inferXBAstrainDynamicsPersistenceRecolonization/inferXBApersistenceRecolonization.R")
rm(list=ls())

# Using the SNP cohorts previously inferred,
# calculate the frequency of each XBA strain.
source("workflow/analysis/inferXBAstrainDynamicsPersistenceRecolonization/inferXBstrainDynamics.R")
rm(list=ls())
