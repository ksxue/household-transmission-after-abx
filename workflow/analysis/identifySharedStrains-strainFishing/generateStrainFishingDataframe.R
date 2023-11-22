# This script aggregates strain from all analyzed species
# and outputs it as a data frame.

library(tidyverse)
library(data.table)
library(foreach)

# Set output directory.
outdir <- "workflow/analysis/identifySharedStrains-strainFishing/out/"

# Import all calculations of the number of private SNPs shared between samples
# for various species in the household cohort.
# Import the dataframe and also annotate it to include the species name.
# This current solution uses the bind_rows function of purrr, via this tutorial:
# https://www.mjandrews.org/blog/readmultifile/
# See also this tutorial on map_dfr: https://dcl-prog.stanford.edu/purrr-extras.html
# See also this post from Claus Wilke, though some of the functions are now outdated:
# https://clauswilke.com/blog/2016/06/13/reading-and-combining-many-tidy-data-files-in-r/
dirs <- Sys.glob("workflow/report/performStrainFishing/*")
files <- sapply(dirs,function(x) paste0(x,"/strainFishing.txt.gz"))
# Check that files exist, since not all directories output strainFishing files.
# These files will not exist if there are no QP samples for that species
# or if there is no list of public SNPs from the HMP.
# (I only have public SNPs for 88 species from the HMP.)
# Note that some species have empty strain fishing output files,
# and they are removed in the snakemake pipeline upstream of this analysis.
filesStrainFishing <- files[file.exists(files)]
# Read in each fixed differences file. Combine the files and annotate by species.
dataStrainFishing <- filesStrainFishing %>%
  map(fread) %>%
  bind_rows(.id="filename") %>%
  separate(filename, into=c("workflow","report","performStrainFishing","species"),
           remove=FALSE, sep="/")
# Keep only columns relevant to downstream analyses.
dataStrainFishing <- dataStrainFishing %>%
  dplyr::select(species, bait, sample, numSitesDetected, numSitesAvailable, meanSNPfreq, medianSNPfreq)


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

# Bind the sample coverage values to the strain fishing dataframe.
dataStrainFishingCov <-
  left_join(dataStrainFishing, dataSampleCoverage %>%
              dplyr::rename(bait=sample, baitcov=cov, baitcovnozeroes=covnozeroes),
            by=c("species","bait"))
dataStrainFishingCov <-
  left_join(dataStrainFishingCov, dataSampleCoverage %>%
              dplyr::rename(samplecov=cov, samplecovnozeroes=covnozeroes),
            by=c("species","sample"))

# Export the dataframe of number of fixed differences.
write.table(dataStrainFishingCov, gzfile(paste0(outdir, "strainFishing-raw.txt.gz")),
            row.names=FALSE, quote=FALSE)
