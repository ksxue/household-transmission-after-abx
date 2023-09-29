library(tidyverse)
library(data.table)
library(foreach)

# Set a global runAll variable that causes all intermediate files to be regenerated.
globalRunAll <- FALSE

# Import list of sequenced samples and condense to sample names.
samplesRaw <- read.table("config/samples-raw.txt", 
                         header=FALSE, stringsAsFactors = FALSE)
colnames(samplesRaw) <- c("sample","raw1","raw2","trimmed1","trimmed2")
# Separate sample names into subject and timepoint.
samplesRaw <- samplesRaw %>%
  dplyr::select(sample) %>% unique() %>%
  separate(sample, into=c("subject","timepoint"), sep="-", remove=FALSE)
samplesRaw$timepoint <- as.numeric(samplesRaw$timepoint)
# Annotate household.
samplesRaw <- samplesRaw %>%
  mutate(hh=substr(sample,1,2))
# Remove the sample with unknown ID (WWW-000).
samplesRaw <- samplesRaw %>%
  filter(!sample=="WWW-000")
# Export list of samples.
write.table(samplesRaw, "workflow/analysis/background/out/samples.txt",
            quote=FALSE, row.names=FALSE)

# Import the list of blacklisted samples and remove them from this list.
sampleBlacklist <- read.table("workflow/analysis/background/out/sampleBlacklist.txt",
                              header=TRUE, stringsAsFactors = FALSE)
samplesRaw <- samplesRaw %>%
  filter(!(sample %in% sampleBlacklist$sample))

# Extract the full list of sequenced samples as well as useful subsets.
samplesAll <- samplesRaw$sample
samplesX <- (samplesRaw %>% filter(substr(sample,1,1)=="X"))$sample
samplesXmain <- (samplesRaw %>% filter(substr(sample,1,1)=="X", timepoint<75))$sample
samplesXfollowup <- (samplesRaw %>% filter(substr(sample,1,1)=="X", timepoint>=75))$sample
samplesY <- (samplesRaw %>% filter(substr(sample,1,1)=="Y"))$sample
samplesZ <- (samplesRaw %>% filter(substr(sample,1,1)=="Z"))$sample
# Extract the list of samples collected at key timepoints.
samplesInitial <- (samplesRaw %>% group_by(subject) %>%
                     filter(timepoint==min(timepoint)))$sample
samplesPreAbx <- sort((samplesRaw %>% group_by(subject) %>%
                    filter(timepoint<30) %>% filter(timepoint==max(timepoint),
                           substr(sample,1,1)=="X"))$sample)
samplesPostAbx <- (samplesRaw %>% group_by(subject) %>%
                    filter(timepoint>35) %>% filter(timepoint==min(timepoint),
                           substr(sample,1,1)=="X"))$sample
samplesFinal <- (samplesRaw %>% group_by(subject) %>%
                   filter(timepoint<75) %>% filter(timepoint==max(timepoint),
                                                   substr(sample,1,1)=="X"))$sample
samplesKeyTimepoints <- 
  c(samplesInitial, samplesPreAbx, samplesPostAbx, samplesFinal)[which(c(samplesInitial, samplesPreAbx, samplesPostAbx, samplesFinal) %in% samplesX)]
samplesKeyTimepointsList <- list(samplesInitial, samplesPreAbx, samplesPostAbx, samplesFinal)
# Extract the list of last samples sequenced for each subject in study arm X.
samplesLastSequencedX <- sort((samplesRaw %>% group_by(subject) %>%
                            filter(substr(sample,1,1)=="X") %>% filter(timepoint==max(timepoint)))$sample)

# Extract the full list of subjects as well as useful subsets.
subjectsAll <- sort(unique(samplesRaw$subject))
subjectsX <- sort(unique((samplesRaw %>% filter(substr(subject,1,1)=="X"))$subject))
# List the antibiotic-taking subjects.
# Based on metagenomic and metabolomic data, annotate XHB as antibiotic-taking
# and XHC as a control subject.
subjectsAbx <- c("XAA","XBA","XCA","XDA","XEA","XFA","XGA","XHB","XIA","XJA",
                 "XKA","XLA","XMA","XNA","XOA","XPA","XQA","XRA","XSA","XTA",
                 "XUA","XVA")
subjectsY <- sort(unique((samplesRaw %>% filter(substr(subject,1,1)=="Y"))$subject))
subjectsZ <- sort(unique((samplesRaw %>% filter(substr(subject,1,1)=="Z"))$subject))
# Extract the list of subjects who participated in follow-up sampling.
subjectsXfollowup <- unique((samplesRaw %>% 
                               filter((timepoint)>75, subject %in% subjectsX))$subject)


# List subjects that appear to switch to alternative states.
subjectsAltStates <- c("XAA","XBA","XDA","XEA","XKA")

# List the additional XBA main study samples that were sequenced later
# to even out the number of timepoints in the study.
samplesXBAextraMain <- c("XBA-041","XBA-042","XBA-045","XBA-046","XBA-047",
                         "XBA-048","XBA-049","XBA-053")
# Subsample the XBA followup timepoints to be approximately monthly,
# like other subjects in the study.
samplesXBAextraFollowup <- 
  c("XBA-351", "XBA-358", "XBA-374", "XBA-380", "XBA-388", "XBA-401", "XBA-408", "XBA-415", 
    #"XBA-422", "XBA-436", "XBA-442", 
    "XBA-455", "XBA-470", "XBA-478", "XBA-482", "XBA-513",
    "XBA-534", "XBA-541", "XBA-557", "XBA-568", "XBA-577", "XBA-589", "XBA-597", "XBA-617",
    "XBA-625", "XBA-650", "XBA-653", "XBA-698", "XBA-706", "XBA-707", "XBA-709", "XBA-710", 
    "XBA-712", "XBA-713", "XBA-716", "XBA-718", "XBA-722", "XBA-727", "XBA-742", "XBA-749", 
    "XBA-754", "XBA-779", "XBA-787", "XBA-794", "XBA-813", "XBA-821", "XBA-828", "XBA-849", 
    "XBA-897", "XBA-898", "XBA-901", "XBA-907")

# Extract the full list of households as well as useful subsets.
hhAll <- sort(unique(samplesRaw$hh))
hhX <- sort(unique((samplesRaw %>% filter(substr(subject,1,1)=="X"))$hh))
hhY <- sort(unique((samplesRaw %>% filter(substr(subject,1,1)=="Y"))$hh))
hhZ <- sort(unique((samplesRaw %>% filter(substr(subject,1,1)=="Z"))$hh))

# Import list of all pairs of cohabiting subjects.
# Note that subjects in households of size >2 will be part of multiple pairs.
cohabitingPairs <- read.table("workflow/analysis/background/cohabitingPairs.txt",
                              header=TRUE, stringsAsFactors = FALSE) %>%
  mutate(pair=paste0(subject1,"-",subject2))
