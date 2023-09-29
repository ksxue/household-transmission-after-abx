# This script takes the species relative abundances
# annotated by antibiotic response and overlays strain sharing info.
# It also summarizes the behaviors for each species in a subject.
library(tidyverse)

# Import metadata about the study.
source("workflow/analysis/background/background.R")
# Set the output directory.
OUTDIR <- "workflow/analysis/annotateSpeciesTrajectoriesSharedStrains/out/"
LIMITOFDETECTION <- 1e-3

# Import annotated species trajectories for all samples in study arm X.
speciesTrajectoriesAll <- 
  read.table(gzfile("workflow/analysis/annotateSpeciesTrajectoriesSharedStrains/out/speciesAbundances-trajectoriesAnnotated-all.txt.gz"),
             header=TRUE, stringsAsFactors = FALSE)
# Import the summary of the annotated species trajectories for all samples in study arm X.
speciesTrajectoriesSummaryAll <- 
  read.table(("workflow/analysis/annotateSpeciesTrajectoriesSharedStrains/out/speciesTrajectoriesSummary-all.txt"),
             header=TRUE, stringsAsFactors = FALSE)


# Import annotated species trajectories for key timepoints and followup samples in study arm X.
speciesTrajectoriesKTfollowup <- 
  read.table(gzfile("workflow/analysis/annotateSpeciesTrajectoriesSharedStrains/out/speciesAbundances-trajectoriesAnnotated-keyTimepointsFollowup.txt.gz"),
             header=TRUE, stringsAsFactors = FALSE)
# Import the summary of the annotated species trajectories for key timepoints and followup samples in study arm X.
speciesTrajectoriesSummaryKTfollowup <- 
  read.table(("workflow/analysis/annotateSpeciesTrajectoriesSharedStrains/out/speciesTrajectoriesSummary-keyTimepointsFollowup.txt"),
             header=TRUE, stringsAsFactors = FALSE)


# Import the species sharing calls within households as a summary per species.
speciesSharingHouseholdsSummary <-
  read.table("workflow/analysis/annotateSpeciesTrajectoriesSharedStrains/out/speciesSharing-hh.txt",
             header=TRUE, stringsAsFactors = FALSE)
# Import the strain sharing calls within households for all timepoints.
strainSharingHouseholdsTimepoint <- 
  read.table("workflow/analysis/compareSharedStrainCalls/out/sameHousehold-sharedStrains-byTimepoint.txt",
             header=TRUE, stringsAsFactors = FALSE)
# Import the strain sharing calls within households as a summary per species.
strainSharingHouseholdsSummary <-
  read.table("workflow/analysis/compareSharedStrainCalls/out/sameHousehold-sharedStrains-subjectSummary.txt",
             header=TRUE, stringsAsFactors = FALSE)


# Import the species sharing calls in the same subject as a summary per species.
speciesSharingSameSubjectSummary <-
  read.table("workflow/analysis/annotateSpeciesTrajectoriesSharedStrains/out/speciesSharing-subject.txt",
             header=TRUE, stringsAsFactors = FALSE)
# Import the strain sharing calls in the same subject as a summary per species.
strainSharingSameSubjectSummary <-
  read.table("workflow/analysis/compareSharedStrainCalls/out/sameSubject-sharedStrains-subjectSummary.txt",
             header=TRUE, stringsAsFactors = FALSE)
# Import the strain sharing calls in the same subject for all timepoints.
strainSharingSameSubjectTimepoint <-
  read.table("workflow/analysis/compareSharedStrainCalls/out/sameSubject-sharedStrains-byTimepoint.txt",
             header=TRUE, stringsAsFactors = FALSE)


# Combine the strain trajectory annotations from all samples and the key timepoints.
speciesTrajectoriesSummaryAnnotated <-
  left_join(speciesTrajectoriesSummaryAll %>% 
              dplyr::rename(speciesDisruptedAllTimepoints=speciesDisrupted,
                            speciesRecoveredAllTimepoints=speciesRecovered,
                            speciesColonizedAllTimepoints=speciesColonized,
                            speciesStatusAbxAllTimepoints=speciesStatusAbx,
                            speciesStatusColonizationAllTimepoints=speciesStatusColonization,
                            timeOfRecoveryAllTimepoints=timeOfRecovery,
                            timeOfColonizationAllTimepoints=timeOfColonization),
            speciesTrajectoriesSummaryKTfollowup %>% 
              dplyr::rename(speciesDisruptedKTfollowup=speciesDisrupted,
                            speciesRecoveredKTfollowup=speciesRecovered,
                            speciesColonizedKTfollowup=speciesColonized,
                            speciesStatusAbxKTfollowup=speciesStatusAbx,
                            speciesStatusColonizationKTfollowup=speciesStatusColonization,
                            timeOfRecoveryKTfollowup=timeOfRecovery,
                            timeOfColonizationKTfollowup=timeOfColonization),
            by=c("subject", "species_id", "kingdom", "phylum", "class", "order",
                 "family", "genus", "species"))

# Annotate the species trajectories based on household strain sharing by subject.
speciesTrajectoriesSummaryAnnotated <- 
  left_join(speciesTrajectoriesSummaryAnnotated, 
            strainSharingHouseholdsSummary %>%
              dplyr::rename(subject=subject1, 
                            strainSharedHhAnyTimepoint=strainSharedAnyTimepoint,
                            strainSharedHhPreAbx=strainSharedPreAbx,
                            strainSharedHhPostAbx=strainSharedPostAbx,
                            strainSharedHhMain=strainSharedMain,
                            strainSharedHhFollowup=strainSharedFollowup,
                            strainSharedHhPostAbxFollowup=strainSharedPostAbxFollowup),
            by=c("species_id", "subject"))
# Annotate the species trajectories based on subject strain sharing.
speciesTrajectoriesSummaryAnnotated <-
  left_join(speciesTrajectoriesSummaryAnnotated,
            strainSharingSameSubjectSummary %>%
              dplyr::rename(strainSharedSubjectInitialPreAbx=strainSharedInitialPreAbx,
                            strainSharedSubjectInitialPostAbx=strainSharedInitialPostAbx,
                            strainSharedSubjectInitialFollowup=strainSharedInitialFollowup,
                            strainSharedSubjectInitialPostAbxFollowup=strainSharedInitialPostAbxFollowup,
                            strainSharedSubjectInitialFinal=strainSharedInitialFinal),
            by=c("species_id", "subject"))
# Annotate the species trajectories based on household species sharing by subject.
speciesTrajectoriesSummaryAnnotated <-
  left_join(speciesTrajectoriesSummaryAnnotated,
            speciesSharingHouseholdsSummary %>% dplyr::select(-maxRelAbundance),
            by=c("species_id", "subject"))
# Annotate the species trajectories based on same-subject species sharing by subject.
speciesTrajectoriesSummaryAnnotated <-
  left_join(speciesTrajectoriesSummaryAnnotated,
            speciesSharingSameSubjectSummary %>% dplyr::select(-maxRelAbundance),
            by=c("species_id", "subject"))

# Annotate species and strain sharing within households for each of the main study periods.
speciesTrajectoriesSummaryAnnotated <-
  speciesTrajectoriesSummaryAnnotated %>%
  mutate(sharingHhAnyTimepoint=
           ifelse(!speciesSharedHhAnyTimepointPresence, "speciesNotShared",
                  ifelse(is.na(strainSharedHhAnyTimepoint), "strainsUnknown",
                         ifelse(strainSharedHhAnyTimepoint, "strainsShared", "strainsNotShared"))),
         sharingHhPreAbx=
           ifelse(!speciesSharedHhPreAbxPresence, "speciesNotShared",
                  ifelse(is.na(strainSharedHhPreAbx), "strainsUnknown",
                         ifelse(strainSharedHhPreAbx, "strainsShared", "strainsNotShared"))),
         sharingHhPostAbx=
           ifelse(!speciesSharedHhPostAbxPresence, "speciesNotShared",
                  ifelse(is.na(strainSharedHhPostAbx), "strainsUnknown",
                         ifelse(strainSharedHhPostAbx, "strainsShared", "strainsNotShared"))),
         sharingHhMain=
           ifelse(!speciesSharedHhMainStudyPresence, "speciesNotShared",
                  ifelse(is.na(strainSharedHhMain), "strainsUnknown",
                         ifelse(strainSharedHhMain, "strainsShared", "strainsNotShared"))),
         sharingHhFollowup=
           ifelse(!speciesSharedHhFollowupPresence, "speciesNotShared",
                  ifelse(is.na(strainSharedHhFollowup), "strainsUnknown",
                         ifelse(strainSharedHhFollowup, "strainsShared", "strainsNotShared"))),
         sharingHhPostAbxFollowup=
           ifelse(!(speciesSharedHhPostAbxPresence | speciesSharedHhFollowupPresence), "speciesNotShared",
                  ifelse(is.na(strainSharedHhPostAbxFollowup), "strainsUnknown",
                         ifelse(strainSharedHhPostAbxFollowup, "strainsShared", "strainsNotShared"))))

# Annotate species and strain sharing within subjects for each of the main study periods.
speciesTrajectoriesSummaryAnnotated <-
  speciesTrajectoriesSummaryAnnotated %>%
  mutate(sharingSubjectPreAbx=
           ifelse(!speciesSubjectPreAbxPresence, "speciesAbsent",
                  ifelse(is.na(strainSharedSubjectInitialPreAbx), "strainsUnknown",
                         ifelse(strainSharedSubjectInitialPreAbx, "strainsShared", "strainsNotShared"))),
         sharingSubjectPostAbx=
           ifelse(!speciesSubjectPostAbxPresence, "speciesAbsent",
                  ifelse(is.na(strainSharedSubjectInitialPostAbx), "strainsUnknown",
                         ifelse(strainSharedSubjectInitialPostAbx, "strainsShared", "strainsNotShared"))),
         sharingSubjectFollowup=
           ifelse(!speciesSubjectFollowupPresence, "speciesAbsent",
                  ifelse(is.na(strainSharedSubjectInitialFollowup), "strainsUnknown",
                         ifelse(strainSharedSubjectInitialFollowup, "strainsShared", "strainsNotShared"))),
         sharingSubjectPostAbxFollowup=
           ifelse(!(speciesSubjectPostAbxPresence | speciesSubjectFollowupPresence), "speciesAbsent",
                  ifelse(is.na(strainSharedSubjectInitialPostAbxFollowup), "strainsUnknown",
                         ifelse(strainSharedSubjectInitialPostAbxFollowup, "strainsShared", "strainsNotShared"))),
         sharingSubjectFinal=
           ifelse(!speciesSubjectFinalSeqPresence, "speciesAbsent",
                  ifelse(is.na(strainSharedSubjectInitialFinal), "strainsUnknown",
                         ifelse(strainSharedSubjectInitialFinal, "strainsShared", "strainsNotShared"))))

# Condense the annotation columns to remove redundant information.
speciesTrajectoriesSummaryAnnotated <- speciesTrajectoriesSummaryAnnotated %>%
  dplyr::select(subject, species_id, kingdom, phylum, class, order, family, genus, species,
                speciesStatusAbxAllTimepoints, speciesStatusColonizationAllTimepoints,
                timeOfRecoveryAllTimepoints, timeOfColonizationAllTimepoints,
                speciesStatusAbxKTfollowup, speciesStatusColonizationKTfollowup,
                timeOfRecoveryKTfollowup, timeOfColonizationKTfollowup, 
                strainTurnover, timeOfStrainTurnover,
                speciesColonizedTransientlyMain, speciesColonizedTransientlyAll,
                starts_with("sharingHh"), starts_with("sharingSubject"))
# Export species trajectory annotations.
write.table(speciesTrajectoriesSummaryAnnotated, 
            paste0(OUTDIR, "speciesTrajectoriesSharingSummary.txt"),
            quote=TRUE, row.names=FALSE)

# Integrate species behavior summaries with full relative abundance data.
speciesTrajectoriesAnnotatedAll <- 
  left_join(speciesTrajectoriesAll %>%
              dplyr::select(sample, subject, hh, timepoint, species_id, 
                            count_reads, coverage, relative_abundance, totalReads),
            speciesTrajectoriesSummaryAnnotated,
            by=c("subject", "species_id"))
# Integrate timepoint-specific information on strain sharing between subjects.
speciesTrajectoriesAnnotatedAll <-
  left_join(speciesTrajectoriesAnnotatedAll,
            strainSharingHouseholdsTimepoint %>%
              dplyr::rename(sample=sample1, sharingHhByTimepoint=strainSharedAnyTimepoint),
            by=c("species_id", "sample"))
# Integrate timepoint-specific information on strain sharing within subjects.
speciesTrajectoriesAnnotatedAll <-
  left_join(speciesTrajectoriesAnnotatedAll,
            strainSharingSameSubjectTimepoint %>%
              dplyr::rename(sharingSubjectTimepointInitial=subjectSharedTimepointInitial,
                            sharingSubjectTimepointFinal=subjectSharedTimepointFinal,
                            sharingSubjectTimepointLastSequenced=subjectSharedTimepointLastSequenced) %>%
              dplyr::select(-strainTurnover, -timeOfStrainTurnover),
            by=c("species_id", "sample"))

# Export species trajectories at all timepoints.
write.table(speciesTrajectoriesAnnotatedAll, 
            gzfile(paste0(OUTDIR, "speciesAbundances-trajectoriesSharingAnnotated-all.txt.gz")),
            quote=TRUE, row.names=FALSE)

