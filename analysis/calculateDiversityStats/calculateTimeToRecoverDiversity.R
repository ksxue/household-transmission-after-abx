library(tidyverse)

# Import background metadata.
source("workflow/analysis/background/background.R")
# Set the output directory.
OUTDIR <- "workflow/analysis/calculateDiversityStats/out/"


# Import the alpha diversity data.
dataAlphaDiversity <- read.table("workflow/analysis/calculateDiversityStats/out/speciesAlpha.txt",
                                 header=TRUE, stringsAsFactors = FALSE) %>%
  mutate(subject=substr(sample,1,3), timepoint=as.integer(substr(sample,5,7)))


# Calculate the median pre-antibiotic diversity for each subject
# and the first timepoint after antibiotics at which it is exceeded.
alphaDiversityRecovery <- dataAlphaDiversity %>%
  filter(subject %in% subjectsAbx) %>%
  group_by(alphaStat, type, subject) %>%
  summarize(medianPreAbxDiversity=median(value[timepoint<29]),
            timeOfRecovery=min(timepoint[timepoint>34 & value>=medianPreAbxDiversity]),
            maxTimepoint=max(timepoint),
            timeOfRecoveryKaplanMeierFormat=ifelse(!is.infinite(timeOfRecovery), as.character(timeOfRecovery),
                     paste0(as.character(maxTimepoint),"+")))

# Export the time to recover the median pre-antibiotic diversity.
write.table(alphaDiversityRecovery, paste0(OUTDIR, "alphaDiversityRecovery.txt"),
            row.names=FALSE, quote=FALSE)
