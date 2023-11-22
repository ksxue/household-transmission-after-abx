# Import data and write a funtion to plot heatmaps of species abundances.

# First import the other supplemental plotting functions.
source("workflow/analysis/plotHelpers.R")

# Plot a heatmap of species abundances in the household of interest.
# Import the annotated species trajectories to plot heatmaps.
# Import the full relative abundances of all annotated species.
dataSpeciesTrajectories <-
  read.table("workflow/analysis/annotateSpeciesTrajectoriesSharedStrains/out/speciesAbundances-trajectoriesSharingAnnotated-all.txt.gz",
             header=TRUE, stringsAsFactors = FALSE)
# Assign each species a short display name.
# Subsample XBA follow-up samples for ease of display.
dataSpeciesTrajectories <- dataSpeciesTrajectories %>%
  ungroup() %>%
  mutate(species_id_short=gsub("_"," ",substr(species_id, 1, nchar(species_id)-6))) %>%
  filter(!(sample %in% samplesXBAextraFollowup), !(sample %in% samplesXBAextraMain))
# Export heatmap summarizing species composition over time.
plotHeatmap <- function(isubject){
  # Focus on a single subject of interest.
  dataSpeciesTrajectoriesSubject <- dataSpeciesTrajectories %>%
    filter(subject==isubject)
  # Annotate the subjects so that species categorized as colonized are not labeled
  # as being disrupted.
  dataSpeciesTrajectoriesSubject <- dataSpeciesTrajectoriesSubject %>%
    mutate(speciesStatusAbxAllTimepoints=
             ifelse(speciesStatusColonizationAllTimepoints=="speciesNotColonized",
                    speciesStatusAbxAllTimepoints, NA))
  # Due to space constraints, retain only species that have median pre-abx relative abundance
  # above the limit of detection or that are annotated as colonized after antibiotics.
  dataSpeciesTrajectoriesSubject <- dataSpeciesTrajectoriesSubject %>%
    group_by(species_id) %>%
    filter((!(is.na(speciesStatusAbxAllTimepoints)) | speciesStatusColonizationAllTimepoints != "speciesNotColonized")) %>%
    filter((median(relative_abundance[timepoint<29])>LIMITOFDETECTION | speciesStatusColonizationAllTimepoints != "speciesNotColonized"))
  # Arrange the species based on response, timing of response, and then alphabetically.
  # Sort the species based on species trajectories.
  dataSpeciesOrdered <- (dataSpeciesTrajectoriesSubject %>%
                           arrange(fct_relevel(speciesStatusAbxAllTimepoints,
                                               c("speciesNotDisrupted", "speciesRecoveredMain",
                                                 "speciesRecoveredFollowup", "speciesNotRecovered")), 
                                   fct_relevel(speciesStatusColonizationAllTimepoints,
                                               c("speciesNotColonized", "pre-abx", "abx", "post-abx", "followup")),
                                   timeOfRecoveryAllTimepoints, 
                                   timeOfColonizationAllTimepoints,
                                   timeOfStrainTurnover) %>%
                           dplyr::select(species_id, species_id_short, speciesStatusAbxAllTimepoints, 
                                         speciesStatusColonizationAllTimepoints) %>% unique() %>%
                           mutate(speciesTrajectory=
                                    ifelse(is.na(speciesStatusAbxAllTimepoints), "colonized",
                                           ifelse(speciesStatusAbxAllTimepoints=="speciesNotDisrupted", "notDisrupted",
                                                  ifelse(speciesStatusAbxAllTimepoints=="speciesNotRecovered", "notRecovered", "recovered")))))
  speciesOrdered <- rev(dataSpeciesOrdered$species_id)
  speciesOrderedShort <- rev(dataSpeciesOrdered$species_id_short)
  speciesColors <- rev(PALETTE.SPECIESTRAJECTORIES[dataSpeciesOrdered$speciesTrajectory])
  # Extract the sequenced timepoints for the subject of interest.
  subjectTimepoints <- sort(unique(dataSpeciesTrajectoriesSubject$timepoint))
  # Extract the relative abundances for samples from the subject of interest
  # and from the initial timepoint in cohabiting partner(s).
  dataSpeciesTrajectoriesHousehold <- dataSpeciesTrajectories %>%
    group_by(subject) %>%
    filter(hh==substr(isubject,1,2),
           species_id %in% dataSpeciesTrajectoriesSubject$species_id) %>%
    ungroup() %>%
    complete(sample, species_id, fill=list(relative_abundance=0)) %>%
    mutate(subject=substr(sample,1,3), timepoint=as.numeric(substr(sample,5,7)))
  dataSpeciesTrajectoriesHousehold <- dataSpeciesTrajectoriesHousehold %>%
    ungroup() %>%
    complete(sample, species_id, fill=list(relative_abundance=0)) %>%
    mutate(subject=substr(sample,1,3), timepoint=as.numeric(substr(sample,5,7))) 
  dataSpeciesTrajectoriesHousehold <- dataSpeciesTrajectoriesHousehold %>%
    mutate(subjectResponse=annotateSubjectResponse(subject)) %>%
    mutate(strainTurnoverAes=ifelse(!is.na(strainTurnover) & strainTurnover & timepoint>=timeOfStrainTurnover,
                                    TRUE, FALSE))
  # Fill in the relative abundances of missing species as 0.
  # Plot a heatmap of the relative abundances of the species present in the community.
  p_heatmap <- dataSpeciesTrajectoriesHousehold %>%
    mutate(timepoint=timepoint-29) %>%
    ggplot() +
    geom_tile(aes(x=factor(timepoint), y=fct_relevel(species_id, speciesOrdered), 
                  fill=log10(relative_abundance))) +
    geom_tile(#data=dataSpeciesTrajectoriesHousehold %>% filter(strainTurnoverAes),
              aes(x=factor(timepoint), y=fct_relevel(species_id, speciesOrdered),
                  color=factor(strainTurnoverAes)),
              alpha=0, linewidth=0.5) +
    geom_point(data=dataSpeciesTrajectoriesHousehold %>%
                 mutate(timepoint=timepoint-29) %>%
                 group_by(subject, species_id) %>%
                 filter(!is.na(sharingHhByTimepoint)) %>%
                 filter(sample %in% c(samplesInitial, samplesLastSequencedX)),
               aes(x=factor(timepoint), y=fct_relevel(species_id, speciesOrdered), 
                   shape=factor(sharingHhByTimepoint)), size=0.8) +
    facet_grid(.~subject, scales="free_x", space="free_x") +
    scale_fill_viridis(option="A", name="Log\nrelative\nabundance") +
    scale_color_manual(name="New\nstrain?", values=c("#FFFFFF00", "#000000")) +
    scale_shape_manual(values=c(1, 16), name="Strain\nshared?") +
    scale_y_discrete(labels=speciesOrderedShort) +
    xlab("Study day") +
    DEFAULTS.THEME_PRINT +
    theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=4.5),
          axis.text.y=element_text(size=4.5, color=speciesColors), 
          strip.text.x=element_text(size=5),
          axis.title.y=element_blank()) +
    guides(fill = guide_colorbar(order=1),
           shape = guide_legend(override.aes = list(size = 1), order=2, rev=TRUE),
           color = guide_legend(override.aes = list(size = 0.5), order=3, rev=TRUE)) +
    theme(legend.key.size = unit(0.5, "lines"))
  p_heatmap
}
