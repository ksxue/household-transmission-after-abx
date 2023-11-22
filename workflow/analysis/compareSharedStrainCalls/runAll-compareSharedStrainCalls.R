# This script consolidates calls of strain sharing from both fixed differences
# and strain sharing into a single set of calls.
# It uses different methods to resolve discrepancies between the two methods.
# The script requires the identifySharedStrains-fixedDiffs and the 
# identifySharedStrains-strainFishing modules to have completed running.

# Consolidate strain sharing calls from fixed differences and strain fishing.
# Annotate discrepancies between the methods and produce a consensus call set
# that uses different methods of resolving the discrepancies.
# Output this set of consolidated strain sharing calls.
source("workflow/analysis/compareSharedStrainCalls/consolidateSharedStrainCalls.R")
rm(list=ls())

# Consolidate strain sharing calls from fixed differences and strain fishing
# for samples from the same subject. This summarizes strain turnover.
# Annotate discrepancies between the methods and produce a consensus call set
# that uses different methods of resolving the discrepancies.
# Output this set of consolidated strain sharing calls.
source("workflow/analysis/compareSharedStrainCalls/consolidateSharedStrainCalls-sameSubject.R")
rm(list=ls())

# Consolidate strain sharing calls from fixed differences and strain fishing
# for samples from the initial timepoint, including both cohabiting and non-cohabiting subjects.
# Annotate discrepancies between the methods and produce a consensus call set
# that uses different methods of resolving the discrepancies.
# Output this set of consolidated strain sharing calls.
source("workflow/analysis/compareSharedStrainCalls/consolidateSharedStrainCalls-sameSubject.R")
rm(list=ls())

# Summarize strain sharing calls from different subjects in the same household 
# by timepoint and by subject.
source("workflow/analysis/compareSharedStrainCalls/summarizeSharedStrainCalls.R")
rm(list=ls())

# Touch a file to verify that the script finished running.
file.create("workflow/analysis/compareSharedStrainCalls/done.txt")
