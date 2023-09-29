library(tidyverse)

# Import metadata about the study.
source("workflow/analysis/background/background.R")
# Set the output directory.
OUTDIR <- "workflow/analysis/calculateTaxonomicEnrichment/out/"
LIMITOFDETECTION <- 1e-3

# Import plot defaults.
source("workflow/analysis/plotDefaults.R")


# Write a function that does enrichment testing ---------------------------


# Write a function to permute species and calculate taxonomic enrichments.
# The function requires a dataframe with group column (usually a taxonomic rank, i.e. family)
# and an annotation column (TRUE/FALSE).
# It also asks for a nickname for the analysis that will be used for outputting
# intermediate data and plots.
calculateEnrichment <- function(dataGroupAnnotation, nickname){
  # Calculate the number of species (rows) with the annotation.
  numSpeciesAnnotated <- nrow(dataGroupAnnotation %>%
                                filter(annotation))
  # Set the seed for this analysis to make it repeatable.
  set.seed(0)
  # Choose a random set of rows with size equal to the number of annotated species.
  # Do this 10000 times.
  randomDraw <- foreach(x=seq(1:10000), .combine="rbind") %do% {
    dataGroupAnnotation %>%
      slice_sample(n=numSpeciesAnnotated) %>%
      dplyr::select(group) %>%
      mutate(replicate=x)
  }
  # Summarize the data from the random draw.
  randomDrawSummary <- randomDraw %>%
    group_by(replicate, group) %>%
    summarize(numSpecies=n()) %>%
    ungroup() %>%
    complete(replicate, group, fill=list(numSpecies=0))
  # Merge the randomized data with the actual data.
  randomDrawSummary <- left_join(randomDrawSummary,
                                 dataGroupAnnotation %>% 
                                   filter(annotation) %>%
                                   group_by(group) %>%
                                   summarize(numSpeciesActual=n()),
                                 by=c("group")) %>%
    mutate(numSpeciesActual=ifelse(is.na(numSpeciesActual),0,numSpeciesActual))
  # Calculate the position of the actual data within the random draw.
  randomDrawPercentiles <- randomDrawSummary %>%
    group_by(group) %>%
    summarize(meanRandom=mean(numSpecies),
              meanActual=mean(numSpeciesActual),
              percentile=sum(numSpecies<=numSpeciesActual)/10000,
              pvalue=ifelse(mean(numSpeciesActual)<mean(numSpecies), 
                            sum(numSpecies<=numSpeciesActual)/10000,
                            sum(numSpecies>=numSpeciesActual)/10000),
              numCounts=n(),
              maxNumSpecies=max(max(numSpecies),meanActual)) %>%
    ungroup() %>%
    mutate(lowerThreshold=0.05/2/n_distinct(group),
           upperThreshold=1-0.05/2/n_distinct(group))
  # Layer the percentile thresholds onto the raw data from the random draws.
  randomDrawSummary <- left_join(randomDrawSummary,
                                 randomDrawPercentiles, by="group")
  
  # Plot the distribution of species in each group in the random draw
  # and the actual number of species.
  p <- randomDrawSummary %>%
    ggplot() +
    geom_boxplot(aes(x=group, y=numSpecies), color="gray50") +
    geom_point(data=randomDrawPercentiles,
               aes(x=group, y=meanActual, color=factor(pvalue<lowerThreshold))) +
    scale_color_manual(values=c("dodgerblue3","firebrick3")) +
    guides(color="none") +
    DEFAULTS.THEME_PRINT +
    xlab("Group") + ylab("Number of species") +
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
  save_plot(paste0(OUTDIR, nickname, "-distribution.png"), p,
            base_width=3, base_height=2.5)
  
  # Export the summary of the random draw and the percentile information.
  write.table(randomDrawPercentiles, paste0(OUTDIR, nickname, "-randomDistributionSummary.txt"),
              row.names=FALSE, quote=FALSE, sep="\t")
  
  # Calculate the likelihood of the actual species distribution.
  dataGroupAnnotationSummary <- dataGroupAnnotation %>% 
    group_by(group, annotation) %>%
    summarize(numSpeciesActual=n()) %>% ungroup() %>%
    complete(group, annotation, fill=list(numSpeciesActual=0))
  dataGroupAnnotationSummary %>%
    mutate(n_ue=numSpeciesActual) %>%
    group_by(group) %>% mutate(n_ue_sumE=sum(n_ue)) %>%
    ungroup() %>% group_by(annotation) %>% mutate(n_ue_sumU=sum(n_ue)) %>%
    ungroup() %>% mutate(n_ue_sumUE=sum(n_ue)) %>%
    mutate(likelihood=n_ue*log((n_ue/n_ue_sumE)/(n_ue_sumU/n_ue_sumUE)))
}


# Perform enrichment testing on species of varying trajectories -----------

# Import the list of species trajectory annotations.
# Annotations include disruption, recovery, and colonization.
speciesTrajectories <- read.table("workflow/analysis/annotateSpeciesTrajectoriesSharedStrains/out/speciesTrajectoriesSharingSummary.txt",
                                  header=TRUE, stringsAsFactors = FALSE)
# Prepare a dataframe to analyze the family-level enrichment among disrupted species
# in the antibiotic-taking subjects.
dataDisrupted <- speciesTrajectories %>%
  filter(!is.na(speciesStatusAbxAllTimepoints),
         subject %in% subjectsAbx) %>%
  dplyr::select(family, speciesStatusAbxAllTimepoints) %>%
  mutate(annotation=(speciesStatusAbxAllTimepoints!="speciesNotDisrupted")) %>%
  dplyr::rename(group=family) %>%
  dplyr::select(group, annotation)
calculateEnrichment(dataDisrupted, "disruptedSpecies")

# Prepare a dataframe to analyze the family-level enrichment among disrupted species
# that do not recover in the antibiotic-taking subjects.
dataDisruptedNoRecovery <- speciesTrajectories %>%
  filter(!is.na(speciesStatusAbxAllTimepoints),
         subject %in% subjectsAbx) %>%
  dplyr::select(family, speciesStatusAbxAllTimepoints) %>%
  mutate(annotation=(speciesStatusAbxAllTimepoints!="speciesNotDisrupted" & 
                       speciesStatusAbxAllTimepoints!="speciesRecoveredMain")) %>%
  dplyr::rename(group=family) %>%
  dplyr::select(group, annotation)
calculateEnrichment(dataDisruptedNoRecovery, "disruptedSpeciesNoRecovery")

# Prepare a dataframe to analyze the family-level enrichment among colonizing species
# and species that experience strain turnovers.
dataColonizationTurnover <- speciesTrajectories %>%
  filter(subject %in% subjectsAbx) %>%
  mutate(annotation=ifelse(speciesStatusColonizationAllTimepoints!="speciesNotColonized" | 
                             (!is.na(strainTurnover) & strainTurnover), TRUE, FALSE)) %>%
  dplyr::select(family, speciesStatusColonizationAllTimepoints, strainTurnover, annotation) %>%
  dplyr::rename(group=family) %>%
  dplyr::select(group, annotation)
calculateEnrichment(dataColonizationTurnover, "colonizingSpeciesStrainTurnover")

# Prepare a dataframe to analyze the family-level enrichment among colonizing species
# and species that experience strain turnovers among antibiotic-taking subjects during the main study.
dataColonizationTurnoverMain <- speciesTrajectories %>%
  filter(subject %in% subjectsAbx) %>%
  mutate(annotation=
           ifelse((speciesStatusColonizationAllTimepoints!="speciesNotColonized" & timeOfColonizationAllTimepoints<75)| 
                     (!is.na(strainTurnover) & strainTurnover & timeOfStrainTurnover<75), TRUE, FALSE)) %>%
  dplyr::select(family, speciesStatusColonizationAllTimepoints, strainTurnover, annotation) %>%
  dplyr::rename(group=family) %>%
  dplyr::select(group, annotation)
calculateEnrichment(dataColonizationTurnoverMain, "colonizingSpeciesStrainTurnover-main")


# Prepare a dataframe to analyze the family-level enrichment among species
# that recover during the main study.
# As a null set, use the full set of species that experienced disruptions.
dataRecovered <- speciesTrajectories %>%
  filter(subject %in% subjectsAbx, 
         !is.na(speciesStatusAbxAllTimepoints),
         speciesStatusAbxAllTimepoints!="speciesNotDisrupted") %>%
  mutate(annotation=ifelse(speciesStatusAbxAllTimepoints=="speciesRecoveredMain",
                           TRUE, FALSE)) %>%
  dplyr::select(family, annotation) %>%
  dplyr::rename(group=family)
calculateEnrichment(dataRecovered, "recoveredSpecies")

# Perform enrichment testing on shared and not shared strains -------------

# Prepare a dataframe to analyze the family-level enrichment among shared strains.
# Import strain sharing annotations at the initial timepoint.
dataStrainSharingInitial <- 
  read.table("workflow/analysis/integrateSpeciesAbundancesSharedStrains/out/annotatedSpeciesAbundances-majorTimepoints.txt.gz",
             header=TRUE, stringsAsFactors = FALSE)

# Import species taxonomies from MIDAS genome database.
speciesTaxons <- fread("workflow/out/midasOutput/database/species_taxonomy.txt",
                       header=TRUE, stringsAsFactors = FALSE, data.table=FALSE)
speciesTaxons <- speciesTaxons %>% 
  dplyr::select(-genome_id,-genome_name,-taxon_id,
                -taxon_lineage_ids,-taxon_lineage_names)
# Shorten some long family names.
speciesTaxons <- speciesTaxons %>%
  mutate(family=gsub("\\..*","",family))
# Add species taxonomy to relative abundances.
dataStrainSharingInitial <- 
  left_join(dataStrainSharingInitial, speciesTaxons,
            by=c("species_id"))

# Extract the set of species in abx-taking subjects that have annotations of strain sharing.
# Set the annotation sample to be the alphabetically first other subject in the household.
# The annotation is true if the strains are shared and false if they are not shared.
# Exclude populations for which strain sharing is not known.
dataStrainSharing <- dataStrainSharingInitial %>%
  mutate(subject=substr(sample,1,3), annotationSubject=substr(annotationSample,1,3)) %>%
  filter(subject %in% subjectsAbx, sample!=annotationSample,
         sample %in% samplesInitial, annotationSample %in% samplesInitial) %>%
  group_by(subject) %>%
  filter(annotationSubject==min(annotationSubject)) %>%
  filter(annotation %in% c("strainsShared", "strainsNotShared")) %>%
  mutate(annotation=(annotation=="strainsShared")) %>% ungroup() %>%
  dplyr::mutate(group=family) %>% dplyr::select(group, annotation)

# Calculate the family-level enrichment of shared strains.
calculateEnrichment(dataStrainSharing, "sharedStrains")
