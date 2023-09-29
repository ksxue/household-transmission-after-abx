# This script loads the strain fishing dataframe and parses out information
# about households and timepoints for downstream analyses.
# It also reads in the list of samples with signatures of cross-contamination
# and removes them from the dataset.

# Import the fixed differences dataframe.
dataStrainFishing <- read.table("workflow/analysis/identifySharedStrains-strainFishing/out/strainFishing-raw.txt.gz",
                                header=TRUE, stringsAsFactors = FALSE)
# Import the list of blacklisted samples.
sampleBlacklist <- (read.table("workflow/analysis/background/out/sampleBlacklist.txt",
                              header=TRUE, stringsAsFactors = FALSE))$sample

# Require that the bait sample (QP) has at least 5x coverage.
dataStrainFishing <- dataStrainFishing %>%
  filter(baitcov>=5)

# Parse the subjects and households from each sample name.
dataStrainFishing <- dataStrainFishing %>%
  mutate(baitSubject=substr(bait,1,3),
         baitHh=substr(bait,1,2),
         baitTimepoint=as.numeric(substr(bait,5,7)),
         sampleSubject=substr(sample,1,3),
         sampleHh=substr(sample,1,2),
         sampleTimepoint=as.numeric(substr(sample,5,7)))
# Annotate comparison as same-subject, same-household, or different-household.
dataStrainFishing <- dataStrainFishing %>%
  mutate(type=ifelse(baitSubject==sampleSubject, "sameSubject",
                     ifelse(baitHh==sampleHh, "sameHousehold", "diffHousehold"))) %>%
  mutate(pctSitesDetected=numSitesDetected/numSitesAvailable)

# Exclude blacklisted samples.
# Remove all calculations in which the bait or sample is on the blacklist.
dataStrainFishing <- dataStrainFishing %>%
  filter(!(bait %in% sampleBlacklist), !(sample %in% sampleBlacklist))
