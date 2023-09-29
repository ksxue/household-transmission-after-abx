# This script aggregates fixed differences data from all sequenced samples
# and outputs it as a data frame.

library(tidyverse)
library(data.table)
library(foreach)

# Set output directory.
outdir <- "workflow/analysis/identifyContaminatedSamples/out/"

# Import all calculations of the number of fixed differences between samples
# for various species in the household cohort.
# Import the dataframe and also annotate it to include the species name.
# This current solution uses the bind_rows function of purrr, via this tutorial:
# https://www.mjandrews.org/blog/readmultifile/
# See also this tutorial on map_dfr: https://dcl-prog.stanford.edu/purrr-extras.html
# See also this post from Claus Wilke, though some of the functions are now outdated:
# https://clauswilke.com/blog/2016/06/13/reading-and-combining-many-tidy-data-files-in-r/
dirs <- Sys.glob("workflow/report/calculateFixedDifferences/*")
files <- sapply(dirs,function(x) paste0(x,"/fixedDiffs_filterCoverage.txt.gz"))
# Check that files exist, since not all directories output fixedDiffs files.
# These files will not exist if there are fewer than 2 samples for that species.
filesFixedDiffs <- files[file.exists(files)]
# Read in each fixed differences file. Combine the files and annotate by species.
dataFixedDiffs <- filesFixedDiffs %>%
  map(fread) %>%
  bind_rows(.id="filename") %>%
  separate(filename, into=c("workflow","report","calculateFixedDifferences","species"),
           remove=FALSE, sep="/") %>%
  dplyr::select(species, sample1, sample2, fixedDiffs, sitesCompared)
# Parse the subjects and households from the information about fixed differences.
dataFixedDiffs <- dataFixedDiffs %>%
  mutate(sample1=substr(sample1, nchar(sample1)-6, nchar(sample1)),
         sample2=substr(sample2, nchar(sample2)-6, nchar(sample2)),
         subject1=substr(sample1,1,3), subject2=substr(sample2,1,3),
         hh1=substr(sample1,1,2), hh2=substr(sample2,1,2),
         timepoint1=as.integer(substr(sample1,5,7)),
         timepoint2=as.integer(substr(sample2,5,7)),
         type=ifelse(subject1==subject2,"sameSubject",
                     ifelse(hh1==hh2,"sameHousehold","diffHousehold")))

# Summarize the number of samples for each species.
dataFixedDiffsSummary <- dataFixedDiffs %>% 
  group_by(species) %>% 
  summarize(numPairs=n(), numSubjects=n_distinct(subject1), 
            numHouseholds=n_distinct(hh1), 
            numWithinHouseholdPairs=sum(type=="sameHousehold")) %>% 
  arrange(desc(numPairs)) %>% 
  mutate(numSamples=(1+sqrt(1+8*numPairs))/2)

# Import all calculations of the median sequencing coverage per sample
# for various species in the household cohort.
# Import the dataframe and also annotate it to include the species name.
filesCoverage <- sapply(dirs,function(x) paste0(x,"/sampleMedianCoverage.txt"))
# Check that files exist.
filesCoverageExists <- filesCoverage[file.exists(filesCoverage)]
# Read in each sample coverage file. Combine the files and annotate by species.
dataSampleCoverage <- filesCoverageExists %>%
  map(fread) %>%
  bind_rows(.id="filename") %>%
  separate(filename, into=c("workflow","report","calculateFixedDifferences","species"),
           remove=FALSE, sep="/") %>%
  dplyr::select(species, sample, cov, covnozeroes) %>%
  mutate(sample=substr(sample, nchar(sample)-6, nchar(sample)))

# Bind the sample coverage values to the fixed differences dataframe.
dataFixedDiffsCov <- 
  left_join(dataFixedDiffs, dataSampleCoverage %>% 
              dplyr::rename(sample1=sample, sample1cov=cov, sample1covnozeroes=covnozeroes),
            by=c("species", "sample1"))
dataFixedDiffsCov <-
  left_join(dataFixedDiffsCov, dataSampleCoverage %>% 
              dplyr::rename(sample2=sample, sample2cov=cov, sample2covnozeroes=covnozeroes),
            by=c("species", "sample2"))

# Export the dataframe of number of fixed differences.
write.table(dataFixedDiffsCov %>%
              dplyr::select(-subject1, -subject2, -hh1, -hh2,
                            -timepoint1, -timepoint2, -type), 
            gzfile(paste0(outdir,"fixedDiffs.txt.gz")),
            row.names=FALSE, quote=FALSE)
# Export the summary of number of samples per species.
write.table(dataFixedDiffsSummary, paste0(outdir,"fixedDiffsSummary.txt"),
            row.names=FALSE, quote=FALSE)
# Export the summary of species coverage from the snps module.
write.table(dataSampleCoverage, paste0(outdir,"sampleCoverage.txt"),
            row.names=FALSE, quote=FALSE)
