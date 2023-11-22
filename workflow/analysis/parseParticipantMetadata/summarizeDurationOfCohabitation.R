library(tidyverse)
library(data.table)

# Set output directory.
OUTDIR <- "workflow/analysis/parseParticipantMetadata/out/"
source("workflow/analysis/background/background.R")

# Import participant metadata.
metadata <- fread("data/householdCohort/ParticipantQuestionnairesCleaned.txt",
                  header=TRUE, data.table=FALSE)

# Extract data on cohabitation.
metadataCohabitation <- metadata %>%
  dplyr::select(`Participant ID`, Housemate1, Housemate2, Housemate3,
                contains("Time living with"), contains("Romantic relationship"))
metadataCohabitationLong <- 
  rbind(metadataCohabitation %>% dplyr::select(`Participant ID`, contains("1")) %>%
          rename_with(function(x){gsub("1","",x)}),
        metadataCohabitation %>% dplyr::select(`Participant ID`, contains("2")) %>%
          rename_with(function(x){gsub("2","",x)}),
        metadataCohabitation %>% dplyr::select(`Participant ID`, contains("3")) %>%
          rename_with(function(x){gsub("3","",x)})) %>%
  filter(!is.na(Housemate)) %>% arrange(`Participant ID`, Housemate)

# Summarize the times of cohabitation and relationship status for each pair.
metadataCohabitationSummary <- metadataCohabitationLong %>%
  dplyr::select(-`Time living with housemate (mo)`) %>%
  mutate(subject1=ifelse(`Participant ID`<Housemate, `Participant ID`, Housemate),
         subject2=ifelse(`Participant ID`<Housemate, Housemate, `Participant ID`)) %>%
  dplyr::select(-`Participant ID`, -Housemate) %>%
  unique() %>%
  filter(subject1!=subject2) %>%
  dplyr::rename(durationCohabitationMonths=`Time living with housemate avg`,
                romanticRelationship=`Romantic relationship housemate`) %>%
  mutate(hh=substr(subject1,1,2)) %>%
  dplyr::select(hh, subject1, subject2, durationCohabitationMonths, romanticRelationship) %>%
  mutate(romanticRelationship=ifelse(romanticRelationship=="Yes",TRUE,FALSE))
# Add back in a row for household XH, which has some missing information.
metadataCohabitationSummary <- 
  rbind(metadataCohabitationSummary,
        c("XH","XHA","XHB",NA,FALSE),
        c("XH","XHB","XHC",NA,TRUE)) %>%
  arrange(subject1, subject2)
  
# Add a variable indicating whether participants have a familial (genetic) relationship
# that is not a romantic relationship.
# In this study, those relationships are parent-child and sibling relationships.
metadataCohabitationSummary <- metadataCohabitationSummary %>%
  mutate(subjectPair=paste0(subject1, "-", subject2)) %>%
  mutate(familialRelationship=
           ifelse(subjectPair %in% c("XAA-XAC", "XAA-XAD", "XHA-XHB", "XHA-XHC",
                                     "XGA-XGC", "XGB-XGC", "XIA-XIC", "XOA-XOB"),
                  TRUE, FALSE)) %>%
  dplyr::select(-subjectPair)

# Export the cleaned summary of cohabitation.
write.table(metadataCohabitationSummary, 
            paste0(OUTDIR, "participantMetadata-cohabitationSummary.txt"),
            row.names=FALSE, quote=FALSE)

# Import the amount of strain sharing between cohabiting subjects at the beginning of the study.
dataRelAbundanceStrainSharingSummaryRaw <-
  read.table(gzfile("workflow/analysis/integrateSpeciesAbundancesSharedStrains/out/annotatedSpeciesAbundances-majorTimepoints-summary.txt"),
             header=TRUE, stringsAsFactors = FALSE)
# Restrict the data to the relevant subjects.
# Also use complete to infer when certain subjects pairs have no strains in common.
dataRelAbundanceStrainSharingSummary <- dataRelAbundanceStrainSharingSummaryRaw %>%
  filter(sample %in% samplesInitial, annotationSample %in% samplesInitial,
         sample!=annotationSample, sample %in% samplesX, annotationSample %in% samplesX,
         substr(sample,1,3)!=substr(annotationSample,1,3)) %>%
  mutate(samplePair=paste0(sample,"_", annotationSample)) %>%
  complete(samplePair, annotation, fill = list(totRelAbundance=0, numSpecies=0)) %>%
  mutate(subject1=substr(samplePair,1,3), subject2=substr(samplePair,9,11),
         subjectPair=ifelse(subject1<subject2, paste0(subject1, "_", subject2), paste0(subject2,"_",subject1))) %>%
  dplyr::select(-subject1, -subject2)
# Parse the data on cohabitation to focus on subjects in study arm X.
metadataCohabitationSummaryX <- metadataCohabitationSummary %>%
  filter(subject1 %in% subjectsX, subject2 %in% subjectsX) %>%
  mutate(subjectPair=paste0(subject1, "_", subject2))
# Combine the data on cohabitation and strain sharing.
dataStrainSharingCohabitation <- 
  left_join(dataRelAbundanceStrainSharingSummary, metadataCohabitationSummaryX,
            by=c("subjectPair"))
dataStrainSharingCohabitation <- dataStrainSharingCohabitation %>%
  mutate(relationshipType=ifelse(romanticRelationship, "romantic",
                          ifelse(familialRelationship, "family", "roommates")))

# Export the combined data on cohabitation and strain sharing.
write.table(dataStrainSharingCohabitation,
            paste0(OUTDIR, "strainSharingInitial-cohabitation.txt"),
            quote=FALSE, row.names=FALSE)
