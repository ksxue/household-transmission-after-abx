# Calculate the JSD from the initial timepoint for each subsequent sample.
# Use the maximum JSD to classify antibiotic-taking subjects into
# minimal, transient, and lasting antibiotic responders.
library(tidyverse)

OUTDIR <- "workflow/analysis/classifySubjectResponses/out/"
# Import study metadata.
source("workflow/analysis/background/background.R")

# Import JSD values calculated from the initial timepoint.
# Import beta diversity data.
dataBeta <- read.table("workflow/analysis/calculateDiversityStats/out/speciesBeta.txt.gz",
                       header=TRUE, stringsAsFactors = FALSE)
# Extract only JSD between the initial timepoint and
# other timepoints from the same subject
dataBetaJSDInitialSameSubject <- dataBeta %>%
  filter(method=="jsd", sample1 %in% samplesInitial,
         substr(sample1,1,3)==substr(sample2,1,3)) %>%
  mutate(subject=substr(sample1,1,3), timepoint=as.numeric(substr(sample2,5,7)))
# Export the JSD values calculated from the initial timepoint for each subject.
write.table(dataBetaJSDInitialSameSubject, paste0(OUTDIR, "speciesJSD-fromInitial-sameSubject.txt"),
            row.names=FALSE, quote=FALSE)

# Extract only JSD between the initial timepoints of subjects
# living in different households.
dataBetaJSDInitialBtwnHh <- dataBeta %>%
  filter(method=="jsd", sample1 %in% samplesInitial,
         substr(sample1,1,2)!=substr(sample2,1,2))
# Export the JSD values calculated between subjects living in different households.
write.table(dataBetaJSDInitialBtwnHh, gzfile(paste0(OUTDIR, "speciesJSD-diffHousehold-initial.txt,gz")),
            row.names=FALSE, quote=FALSE)
# Calculate the 95th percentile range of JSD values between unrelated subjects.
JSDthresholds <- quantile(dataBetaJSDInitialBtwnHh$value,
                          c(0.025,0.25,0.75,0.975))

# Calculate the maximum JSD from the initial timepoint for each subject
# and the final JSD at the end of the main study.
dataMaxFinalJSD <- dataBetaJSDInitialSameSubject %>%
  filter(sample2 %in% samplesXmain) %>%
  group_by(subject) %>%
  summarize(maxJSD=max(value), finalJSD=value[timepoint==max(timepoint)])
# Classify each subject response based on the maximum JSD.
dataSubjectResponses <- dataMaxFinalJSD %>%
  mutate(subjectResponse=
           ifelse(!(subject %in% subjectsAbx), "control",
           ifelse(maxJSD<0.4, "minimal response",
           ifelse(finalJSD<0.4, "transient response", "lasting response"))))
# Export the subject classifications.
write.table(dataSubjectResponses %>% dplyr::select(subject, subjectResponse), 
            paste0(OUTDIR, "subjectResponses.txt"),
            row.names=FALSE, quote=FALSE, sep="\t")
# Export the maximum and final JSD.
write.table(dataMaxFinalJSD, paste0(OUTDIR, "subject-maxFinalJSD.txt"),
            row.names=FALSE, quote=FALSE, sep="\t")
