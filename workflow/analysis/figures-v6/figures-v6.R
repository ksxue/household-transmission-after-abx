

# Figure 1 ----------------------------------------------------------------

# Import plot defaults and helper functions.
source("workflow/analysis/plotHelpers.R")
# Set output directory.
OUTDIR <- "workflow/analysis/figures-v6/out/"

fig1summarywidth <- 1.25
fig1height <- 1.25
fig1relativeabundancewidth <- 1.8
fig1trajectorywidth <- 3.75

# Calculate the median number of follow-up samples for each subject
# and the median duration of follow-up sampling.
# Do not include subjects who collected zero follow-up samples.
# median: 3 samples; median duration: 364 days
# num subjects > 2 years: 6
samplesRaw %>%
  filter(!(sample %in% sampleBlacklist), timepoint>75, subject %in% subjectsX) %>%
  group_by(subject) %>% summarize(numSamples=n(), followupDuration=max(timepoint)-64) %>%
  ungroup() %>% summarize(medianNumSamples=median(numSamples), medianDurationDays=median(followupDuration),
                          medianDurationYears=medianDurationDays/365, medianDurationMonths=medianDurationYears*12,
                          numSubjectsMoreThanTwoYears=sum(followupDuration>2*365), numSubjects=n(),
                          totalSamplingYrs=sum(followupDuration+64)/365,
                          totalSamplingYrsNoMainStudy=sum(followupDuration)/365)
# Compare the sampling in subjects with lasting antibiotic responses compared to
# other subjects.
samplesRaw %>%
  filter(!(sample %in% sampleBlacklist), timepoint>75, subject %in% subjectsX) %>%
  group_by(subject) %>%
  summarize(numSamples=n(), followupDuration=max(timepoint)-64) %>% ungroup() %>%
  mutate(subjectResponse=annotateSubjectResponse(subject)) %>%
  group_by(subjectResponse) %>%
  summarize(medianNumSamples=median(numSamples), medianDurationDays=median(followupDuration))

# Plot the relative abundance of shared species and strains
# at the initial timepoint for subjects in study arm X.
# Import the relative abundance data annotated by strain sharing.
dataRelAbundanceStrainSharingRaw <-
  read.table(gzfile("workflow/analysis/integrateSpeciesAbundancesSharedStrains/out/annotatedSpeciesAbundances-majorTimepoints.txt.gz"),
             header=TRUE, stringsAsFactors = FALSE)
# Import the summary of the relative abundance data.
dataRelAbundanceStrainSharingSummaryRaw <-
  read.table(gzfile("workflow/analysis/integrateSpeciesAbundancesSharedStrains/out/annotatedSpeciesAbundances-majorTimepoints-summary.txt"),
             header=TRUE, stringsAsFactors = FALSE)
# Restrict the data to the relevant subjects.
dataRelAbundanceStrainSharing <- dataRelAbundanceStrainSharingRaw %>%
  filter(sample %in% samplesInitial, annotationSample %in% samplesInitial,
         sample!=annotationSample, sample %in% samplesX, annotationSample %in% samplesX,
         substr(sample,1,3)!=substr(annotationSample,1,3)) %>%
  mutate(samplePair=paste0(sample,"_",annotationSample))
# Also use complete to infer when certain subjects pairs have no strains in common.
dataRelAbundanceStrainSharingSummary <- dataRelAbundanceStrainSharingSummaryRaw %>%
  filter(sample %in% samplesInitial, annotationSample %in% samplesInitial,
         sample!=annotationSample, sample %in% samplesX, annotationSample %in% samplesX,
         substr(sample,1,3)!=substr(annotationSample,1,3)) %>%
  mutate(samplePair=paste0(sample,"_", annotationSample)) %>%
  complete(samplePair, annotation, fill = list(totRelAbundance=0, numSpecies=0))

# Calculate the maximum relative abundance of shared strains in any cohabiting pair.
# Most recent max value: 66.3%, XGB and XGC
# Most recent values: min 0, max 66%, median 7.6%
dataRelAbundanceStrainSharingSummary %>%
  filter(annotation=="strainsShared") %>%
  summarize(medianAbundance=median(totRelAbundance),
            minAbundance=min(totRelAbundance), maxAbundance=max(totRelAbundance))

# Calculate the maximum number of shared strains in any cohabiting pair.
# Most recent values: min 0, max 12, median 3
dataRelAbundanceStrainSharingSummary %>%
  filter(annotation=="strainsShared") %>%
  summarize(medianAbundance=median(numSpecies),
            minAbundance=min(numSpecies), maxAbundance=max(numSpecies))

# Calculate the minimum and maximum relative abundance of non-shared strains in any cohabiting pair.
# Most recent values: min 1.6%, median 36%, max 71%
dataRelAbundanceStrainSharingSummary %>%
  filter(annotation=="strainsNotShared") %>%
  summarize(medianAbundance=median(totRelAbundance),
            minAbundance=min(totRelAbundance), maxAbundance=max(totRelAbundance))

# Calculate the number of households in which we find a cohabiting pair
# that shares zero strains (note that other pairs in the household may share strains).
# Most recent value: 5 households, 22.7% (5/22)
dataRelAbundanceStrainSharingSummary %>%
  filter(annotation=="strainsShared", numSpecies==0) %>%
  mutate(hh=substr(samplePair,1,2)) %>% pull(hh) %>% unique() %>% length()
# Calculate the number of cohabiting pairs that have zero strain sharing.
# Most recent value: 7 cohabiting pairs, 23% (7/30)
dataRelAbundanceStrainSharingSummary %>%
  filter(annotation=="strainsShared", numSpecies==0) %>%
  mutate(subject1=substr(samplePair,1,3), subject2=substr(samplePair,9,11),
         cohabitingPair=ifelse(subject1<subject2,
                               paste0(subject1,"-",subject2), paste0(subject2,"-",subject1))) %>%
  pull(cohabitingPair) %>% unique() %>% length()

# Calculate the median non-zero abundance of shared strains
# Most recent value: min 0.8%, median 23%, max 66%
dataRelAbundanceStrainSharingSummary %>%
  filter(annotation=="strainsShared", numSpecies!=0) %>%
  summarize(medianNonzeroSharing=median(totRelAbundance),
            minNonzeroSharing=min(totRelAbundance),
            maxNonzeroSharing=max(totRelAbundance))
# Calculate the maximum number of shared strains in cohabiting pairs that share strains.
# Most recent values: min 0, max 12, median 4
dataRelAbundanceStrainSharingSummary %>%
  filter(annotation=="strainsShared", numSpecies!=0) %>%
  summarize(medianNumSpecies=median(numSpecies),
            minNumSpecies=min(numSpecies), maxNumSpecies=max(numSpecies))

# Plot the distribution of strains that are shared and not shared
# at the beginning of the study.
p1_sharingSummaryHistogram <- dataRelAbundanceStrainSharingSummary %>%
  filter(annotation %in% c("strainsShared","strainsNotShared")) %>%
  mutate(displayAnnotation=
           ifelse(annotation=="strainsShared","species shared,\nstrains shared",
                  "species shared,\nstrains not shared")) %>%
  ggplot() +
  geom_histogram(aes(x=totRelAbundance, fill=factor(annotation)), binwidth=0.05) +
  facet_wrap(~fct_rev(displayAnnotation)) +
  scale_fill_manual(values=sharingAnnotationPalette,
                    name="Strain\nshared?", labels=c(FALSE, TRUE)) +
  xlab("Total abundance") + ylab("Number of\ncohabiting pairs") +
  guides(fill = "none") + xlim(-0.05,1) +
  DEFAULTS.THEME_PRINT + theme(strip.text.x=element_text(size=6))
savePNGPDF(paste0(OUTDIR, "1-strainSharingSummaryHistogram"), p1_sharingSummaryHistogram,
           fig1height, 2.5)

# Plot the relationship between time of cohabitation and proportion of strains shared.
# Import the data on strain sharing and the duration of cohabitation.
dataStrainSharingCohabitation <- read.table("workflow/analysis/parseParticipantMetadata/out/strainSharingInitial-cohabitation.txt",
                                            header=TRUE, stringsAsFactors = FALSE) %>%
  filter(annotation=="strainsShared")
# Calculate the number of subject pairs who shared >30% of their gut microbiomes
# despite living together for less than one year.
# 3 subject pairs share >30% of their gut microbiomes despite living together
# for less than one year.
# Do not include household XI, since familial relationships could confound this relationship.
dataStrainSharingCohabitation %>%
  filter(durationCohabitationMonths<12, totRelAbundance>0.3) %>%
  dplyr::select(subjectPair) %>% unique() %>% nrow()
# Calculate the Pearson correlation between length of cohabitation (log 10)
# and the amount of strain sharing.
# Do not use Spearman, since each pair of subjects has the same length
# of cohabitation, which creates problems with calculating rank.
# Pearson r=0.3, p=0.02075
cohabitationStrainSharing_pearson <-
  cor.test((dataStrainSharingCohabitation %>%
              mutate(logDurationCohabitationYears=log10(durationCohabitationMonths/12)))$logDurationCohabitationYears,
           (dataStrainSharingCohabitation)$totRelAbundance,
           method="pearson")
# Test if the correlation holds without the log10.
# Pearson r=0.21, p=0.12 (not significant)
cohabitationStrainSharing_pearson_noLog <-
  cor.test((dataStrainSharingCohabitation)$durationCohabitationMonths,
           (dataStrainSharingCohabitation)$totRelAbundance,
           method="pearson")
# Plot the relationship between duration of cohabitation and the percent
# of the microbiome shared.
p1_strainSharingCohabitation <- dataStrainSharingCohabitation %>%
  ggplot() +
  geom_point(aes(x=durationCohabitationMonths, y=totRelAbundance,
                 color=factor(relationshipType)), alpha=0.8, size=1) +
  geom_text(aes(x=0.3, y=0.9,
                label=paste0("Pearson r= ", round(cohabitationStrainSharing_pearson$estimate, digits=2),
                             "\np= ", round(cohabitationStrainSharing_pearson$p.value, digits=2))),
            size=2, hjust=0) +
  scale_color_manual(values=DEFAULTS.PALETTE.TOL.COLORBLINDSAFE, name="Relationship") +
  scale_x_log10() + ylim(0,1) +
  xlab("Length of cohabitation\n(months)") +
  ylab("Total abundance,\nshared strains") +
  DEFAULTS.THEME_PRINT +
  guides(color = guide_legend(override.aes = list(size = 0.75))) +
  theme(legend.key.size = unit(0.5, "lines")) +
  theme(legend.box.margin=margin(0,0,0,-10))
savePNGPDF(paste0(OUTDIR, "1-strainSharingCohabitation"), p1_strainSharingCohabitation,
           fig1height+0.1, 2.2)

# Calculate the mean amount of strain sharing and length of cohabitation
# for each household.
# Pearson p=0.41
dataStrainSharingCohabitationByHouseholdMean <- dataStrainSharingCohabitation %>%
  filter(!is.na(durationCohabitationMonths)) %>%
  group_by(hh) %>%
  summarize(meanProportionSharing=mean(totRelAbundance),
         meanDurationCohabitation=mean(durationCohabitationMonths))
cor.test(dataStrainSharingCohabitationByHouseholdMean$meanProportionSharing,
         log10(dataStrainSharingCohabitationByHouseholdMean$meanDurationCohabitation),
         method="pearson")
# For pairs, calculate the mean amount of strain sharing and length of cohabitation
# for each household. For trios, take the maximum length of cohabitation.
# Pearson p=0.68
dataStrainSharingCohabitationByHouseholdMax <- dataStrainSharingCohabitation %>%
  filter(!is.na(durationCohabitationMonths)) %>%
  group_by(hh) %>%
  filter(n()==2 | durationCohabitationMonths==max(durationCohabitationMonths)) %>%
  summarize(meanProportionSharing=mean(totRelAbundance),
            meanDurationCohabitation=mean(durationCohabitationMonths),
            relationshipType=unique(relationshipType))
cohabitationStrainSharing_pearson_byHouseholdMax <- 
  cor.test(dataStrainSharingCohabitationByHouseholdMax$meanProportionSharing,
         log10(dataStrainSharingCohabitationByHouseholdMax$meanDurationCohabitation),
         method="pearson")
# For pairs, calculate the mean amount of strain sharing and length of cohabitation
# for each household. For trios, take the minimum length of cohabitation.
# Pearson p=0.68
dataStrainSharingCohabitationByHouseholdMin <- dataStrainSharingCohabitation %>%
  filter(!is.na(durationCohabitationMonths)) %>%
  group_by(hh) %>%
  filter(n()==2 | durationCohabitationMonths==min(durationCohabitationMonths)) %>%
  summarize(meanProportionSharing=mean(totRelAbundance),
            meanDurationCohabitation=mean(durationCohabitationMonths))
cor.test(dataStrainSharingCohabitationByHouseholdMin$meanProportionSharing,
         log10(dataStrainSharingCohabitationByHouseholdMin$meanDurationCohabitation),
         method="pearson")

# Plot the relationship between duration of cohabitation and the percent
# of the microbiome shared.
p1_strainSharingCohabitationByHousehold <- dataStrainSharingCohabitationByHouseholdMax %>%
  ggplot() +
  geom_point(aes(x=meanDurationCohabitation, y=meanProportionSharing,
                 color=factor(relationshipType)), alpha=0.8, size=1) +
  annotate(geom="text", x=0.3, y=0.9,
                label=paste0("Pearson r= ", round(cohabitationStrainSharing_pearson_byHouseholdMax$estimate, digits=2),
                             "\np= ", round(cohabitationStrainSharing_pearson_byHouseholdMax$p.value, digits=2)),
            size=2, hjust=0) +
  scale_color_manual(values=DEFAULTS.PALETTE.TOL.COLORBLINDSAFE, name="Relationship") +
  scale_x_log10() + ylim(0,1) +
  xlab("Length of cohabitation\n(months)") +
  ylab("Total abundance,\nshared strains") +
  DEFAULTS.THEME_PRINT +
  guides(color = guide_legend(override.aes = list(size = 0.75))) +
  theme(legend.key.size = unit(0.5, "lines")) +
  theme(legend.box.margin=margin(0,0,0,-10))
savePNGPDF(paste0(OUTDIR, "1-strainSharingCohabitation-byHousehold"), p1_strainSharingCohabitationByHousehold,
           fig1height+0.1, 2.2)

# Import species abundances and defaults for relative abundance plots.
source("workflow/analysis/generateSpeciesAbundances/loadSpeciesAbundancePlots.R")

# Plot community composition during the main study for each subject.
# Plot a control subject.
savePNGPDF(paste0(OUTDIR, "1-communityComposition-control-XIC"),
           plotSpeciesAbundanceDaysFromAbxStart(speciesAbundances %>%
                                  filter(subject=="XIC", sample %in% samplesXmain) %>%
                                  group_by(species_id) %>% filter(max(relative_abundance)>LIMITOFDETECTION) %>%
                                  mutate(subject="")) +
             DEFAULTS.THEME_PRINT + ylab("Relative\nabundance"), fig1height, fig1relativeabundancewidth)
# Plot a subject with a minimal antibiotic response.
savePNGPDF(paste0(OUTDIR, "1-communityComposition-minimal-XVA"),
           plotSpeciesAbundanceDaysFromAbxStart(speciesAbundances %>%
                                  filter(subject=="XVA", sample %in% samplesXmain) %>%
                                  group_by(species_id) %>% filter(max(relative_abundance)>LIMITOFDETECTION) %>%
                                  mutate(subject="")) +
             DEFAULTS.THEME_PRINT + ylab("Relative\nabundance"), fig1height, fig1relativeabundancewidth)
# Plot a subject with a transient antibiotic response.
savePNGPDF(paste0(OUTDIR, "1-communityComposition-transient-XSA"),
           plotSpeciesAbundanceDaysFromAbxStart(speciesAbundances %>%
                                  filter(subject=="XSA", sample %in% samplesXmain) %>%
                                  group_by(species_id) %>% filter(max(relative_abundance)>LIMITOFDETECTION) %>%
                                  mutate(subject="")) +
             DEFAULTS.THEME_PRINT + ylab("Relative\nabundance"), fig1height, fig1relativeabundancewidth)
# Plot a subject with a lasting antibiotic response.
savePNGPDF(paste0(OUTDIR, "1-communityComposition-lasting-XDA"),
           plotSpeciesAbundanceDaysFromAbxStart(speciesAbundances %>%
                                  filter(subject=="XDA", sample %in% samplesXmain) %>%
                                  group_by(species_id) %>% filter(max(relative_abundance)>LIMITOFDETECTION) %>%
                                  mutate(subject="")) +
             DEFAULTS.THEME_PRINT + ylab("Relative\nabundance"), fig1height, fig1relativeabundancewidth)

# Export the common families palette.
# List only the top 15 families.
topFamilies <- (speciesAbundances %>%
                  filter(sample %in% samplesInitial) %>%
                  group_by(family) %>% summarize(totAbundance=sum(relative_abundance)) %>%
                  top_n(15, totAbundance) %>%
                  ungroup() %>% mutate(family=ifelse(family=="","Unknown",family)))$family
topFamiliesPalette <- commonFamiliesPalette[topFamilies]
topFamiliesdataframe <- data.frame(family=topFamilies, color=topFamiliesPalette)
topFamiliesdataframe$family <- factor(topFamilies,
                                      levels=c("Unknown", topFamilies[-which(topFamilies=="Unknown")]))
p1_familiesPalette <- topFamiliesdataframe %>%
  ggplot() +
  geom_tile(aes(x=family, y=0, fill=factor(family))) +
  scale_fill_manual(name="Family", values=topFamiliesPalette) +
  DEFAULTS.THEME_PRINT +
  guides(fill = guide_legend(override.aes = list(size = 0.5))) +
  theme(legend.key.size = unit(0.5, "lines"))
savePNGPDF(paste0(OUTDIR, "1-relativeAbundance-palette"), get_legend(p1_familiesPalette),
           1.5*fig1height, 1)
p1_familiesPalette <- topFamiliesdataframe %>%
  ggplot() +
  geom_tile(aes(x=family, y=0, fill=factor(family))) +
  scale_fill_manual(name="Family", values=topFamiliesPalette) +
  DEFAULTS.THEME_PRINT +
  guides(fill = guide_legend(override.aes = list(size = 0.5), ncol=2)) +
  theme(legend.key.size = unit(0.5, "lines"))
savePNGPDF(paste0(OUTDIR, "1-relativeAbundance-palette-2col"),
           get_legend(p1_familiesPalette +
                        guides(fill = guide_legend(override.aes = list(size = 0.5), ncol=2))),
           1.25, 2)
savePNGPDF(paste0(OUTDIR, "1-relativeAbundance-palette-3col"),
           get_legend(p1_familiesPalette +
                        guides(fill = guide_legend(override.aes = list(size = 0.5), ncol=3))),
           1.25, 3)


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
# Extract only JSD between the initial timepoints of subjects
# living in different households.
dataBetaJSDInitialBtwnHh <- dataBeta %>%
  filter(method=="jsd", sample1 %in% samplesInitial,
         substr(sample1,1,2)!=substr(sample2,1,2))
JSDthresholds <- quantile(dataBetaJSDInitialBtwnHh$value,
                          c(0.025,0.25,0.75,0.975))

# Plot the JSD trajectories of individual subjects.
plotJSDtrajectoryDaysFromAbxStart <- function(x){
  dataBetaJSDInitialSameSubject %>%
    filter(subject==x) %>%
    mutate(timepoint=timepoint-29) %>%
    mutate(subjectResponse=annotateSubjectResponse(subject), display="") %>%
    ggplot() +
    geom_rect(aes(xmin=-29, xmax=35, ymin=JSDthresholds[1], ymax=JSDthresholds[4]),
              fill="gray70") +
    geom_rect(aes(xmin=0, xmax=5, ymin=0, ymax=1),
              fill=ifelse(x %in% subjectsAbx, abxColor, NA), alpha=0.8) +
    geom_line(aes(x=timepoint, y=value, group=subject, color=factor(subjectResponse))) +
    facet_wrap(~display) +
    scale_color_manual(values=PALETTE.ABXRESPONSE) +
    scale_x_continuous(breaks=timepointScaleBreaks, limits=c(-29,35)) +
    ylim(0,1) +
    xlab("Study day") + ylab("JSD from\nstudy start") + guides(color="none") +
    DEFAULTS.THEME_PRINT
}
savePNGPDF(paste0(OUTDIR, "1-JSD-control-XIC"), plotJSDtrajectoryDaysFromAbxStart("XIC"), fig1height, fig1height+0.25)
savePNGPDF(paste0(OUTDIR, "1-JSD-minimal-XVA"), plotJSDtrajectoryDaysFromAbxStart("XVA"), fig1height, fig1height+0.25)
savePNGPDF(paste0(OUTDIR, "1-JSD-transient-XSA"), plotJSDtrajectoryDaysFromAbxStart("XSA"), fig1height, fig1height+0.25)
savePNGPDF(paste0(OUTDIR, "1-JSD-lasting-XDA"), plotJSDtrajectoryDaysFromAbxStart("XDA"), fig1height, fig1height+0.25)


# Plot the maximum and final JSD for each subject.
p1_JSDmaxfinal <- dataBetaJSDInitialSameSubject %>%
  filter(sample2 %in% samplesXmain) %>%
  group_by(subject) %>%
  summarize(maxJSD=max(value),
            finalJSD=value[timepoint==max(timepoint)]) %>%
  mutate(subjectResponse=annotateSubjectResponse(subject),
         subjectResponse=gsub(" response","",subjectResponse)) %>%
  ggplot() +
  #geom_vline(aes(xintercept=0.4), linetype="dashed") +
  geom_hline(aes(yintercept=0.4), linetype="dashed") +
  #geom_abline(intercept=0, slope=1, linetype="dashed") +
  geom_point(aes(x=maxJSD, y=finalJSD,
                 color=fct_relevel(subjectResponse, ABXRESPONSESHORT)), alpha=0.8, size=1) +
  scale_color_manual(values=PALETTE.ABXRESPONSESHORT, name="Subject\nresponse") +
  xlim(0,1) + ylim(0,1) +
  guides(color = guide_legend(override.aes = list(size = 0.75))) +
  theme(legend.key.size = unit(0.5, "lines")) +
  xlab("Maximum JSD\nfrom study start") + ylab("JSD to study start,\n30 days post-abx") +
  DEFAULTS.THEME_PRINT +
  theme(legend.box.margin=margin(0,0,0,-10))
savePNGPDF(paste0(OUTDIR, "1-JSD-max-final"), p1_JSDmaxfinal,
           fig1height+0.1, 2.15)


# Fig 2 - Strain and species dynamics after antibiotic perturbation -------

# Import plot defaults and helper functions.
source("workflow/analysis/plotHelpers-heatmaps.R")
# Set output directory.
OUTDIR <- "workflow/analysis/figures-v6/out/"

# Plot heatmap for household XD, focusing on subject XDA.
p2_heatmapXDA <- plotHeatmap("XDA") + guides(fill="none", color="none", shape="none")
savePNGPDF(paste0(OUTDIR, "2-heatmap-XD"),p2_heatmapXDA, 4.5, 5)
# Export the legend for the heatmaps.
savePNGPDF(paste0(OUTDIR, "2-heatmap-legend"), get_legend(plotHeatmap("XDA")), 2.5, 0.75)

# Plot heatmap for household XA, focusing on subject XAA.
savePNGPDF(paste0(OUTDIR, "S7-heatmap-XA"),
           plotHeatmap("XAA") + guides(fill="none",color="none",shape="none"), 5, 5)
savePNGPDF(paste0(OUTDIR, "S7-heatmap-XA-legend"), get_legend(plotHeatmap("XAA")), 2.5, 0.75)

# Plot heatmap for household XB, focusing on subject XBA.
savePNGPDF(paste0(OUTDIR, "S8-heatmap-XB"),
           plotHeatmap("XBA") + guides(fill="none",color="none",shape="none"), 4.5, 6.25)
savePNGPDF(paste0(OUTDIR, "S8-heatmap-XB-legend-fill"),
           get_legend(plotHeatmap("XBA") + guides(color="none",shape="none")), 1, 0.75)
savePNGPDF(paste0(OUTDIR, "S8-heatmap-XB-legend-color"),
           get_legend(plotHeatmap("XBA") + guides(fill="none",shape="none")), 1, 0.75)
savePNGPDF(paste0(OUTDIR, "S8-heatmap-XB-legend-shape"),
           get_legend(plotHeatmap("XBA") + guides(fill="none",color="none")), 1, 0.75)


# Plot heatmap for household XE, focusing on subject XEA.
savePNGPDF(paste0(OUTDIR, "S9-heatmap-XE"),
           plotHeatmap("XEA") + guides(fill="none",color="none",shape="none"), 5, 4)
savePNGPDF(paste0(OUTDIR, "S9-heatmap-XE-legend"), get_legend(plotHeatmap("XEA")), 2.5, 0.75)

# Plot heatmap for household XK, focusing on subject XKA.
savePNGPDF(paste0(OUTDIR, "S10-heatmap-XK"),
           plotHeatmap("XKA") + guides(fill="none",color="none",shape="none"), 4, 5)
savePNGPDF(paste0(OUTDIR, "S10-heatmap-XK-legend"), get_legend(plotHeatmap("XKA")), 2.5, 0.75)


p1_familiesPalette <- topFamiliesdataframe %>%
  ggplot() +
  geom_tile(aes(x=family, y=0, fill=factor(family))) +
  scale_fill_manual(name="Family", values=topFamiliesPalette) +
  DEFAULTS.THEME_PRINT +
  guides(fill = guide_legend(override.aes = list(size = 0.5))) +
  theme(legend.key.size = unit(0.5, "lines"))
savePNGPDF(paste0(OUTDIR, "1-relativeAbundance-palette"), get_legend(p1_familiesPalette),
           1.5*fig1height, 1)

# Export the trajectories of species of interest.
plotSpeciesTrajectory <- function(ihh, ispecies_id, height, width){
  # Extract data for species of interest.
  dataSpecies <- dataSpeciesTrajectories %>%
    filter(hh==ihh, species_id==ispecies_id) %>%
    mutate(relative_abundance=ifelse(relative_abundance<1e-5, 1e-5, relative_abundance)) %>%
    mutate(timePeriod=ifelse(sample %in% samplesXmain, "main", "followup")) %>%
    mutate(subjectResponse=gsub(" ","\n", annotateSubjectResponse(subject)))
  # Set a color for the abx subject based on the color of subject responses.
  abxSubject <- unique((dataSpecies %>% filter(subject %in% subjectsAbx))$subject)
  strainSharingShapes <- c(1,16)
  names(strainSharingShapes) <- c(FALSE,TRUE)
  abxStartIndex <- which(sort(unique(dataSpecies$timepoint))==29)
  abxEndIndex <- which(sort(unique(dataSpecies$timepoint))==34)
  followupStartIndex <- min(which(sort(unique(dataSpecies$timepoint))>75))
  followupEndIndex <- max(which(sort(unique(dataSpecies$timepoint))>75))
  p2_species <- dataSpecies %>%
    arrange((subjectResponse)) %>%
    ggplot() +
    scale_y_continuous(trans='log10', limits=c(1e-5,1),
                       breaks=trans_breaks('log10', function(x) 10^x),
                       labels=trans_format('log10', math_format(10^.x))) +
    scale_x_discrete() +
    geom_hline(yintercept=1e-5, linetype="dashed") +
    annotate(geom="rect", xmin=abxStartIndex, xmax=abxEndIndex, ymin=10^-5, ymax=1,
              fill=abxColor, alpha=0.5) +
    annotate(geom="rect", xmin=followupStartIndex, xmax=followupEndIndex, ymin=10^-5, ymax=1,
              fill="grey80") +
    geom_line(aes(x=factor(timepoint-29), y=relative_abundance, group=subject,
                  color=factor(subjectResponse))) +
    geom_point(data=dataSpecies %>% filter(!is.na(sharingHhByTimepoint)) %>%
                 arrange(subjectResponse),
               aes(x=factor(timepoint-29), y=relative_abundance,
                   color=factor(subjectResponse),
                   shape=factor(sharingHhByTimepoint)), alpha=0.8, size=1) +
    scale_color_manual(values=PALETTE.ABXRESPONSELINEBREAK, name="Subject") +
    scale_shape_manual(values=strainSharingShapes, name="Strain\nshared?") +
    xlab("Study day") + ylab("Relative\nabundance") +
    DEFAULTS.THEME_PRINT +
    ggtitle(gsub("_"," ",substr(ispecies_id, 1, nchar(ispecies_id)-6))) +
    theme(plot.title=element_text(size=5, margin=margin(0,0,1,0)),
          strip.text.x=element_blank()) +
    theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=3),
          axis.text.y=element_text(size=4),
          axis.title=element_text(size=5)) +
    guides(shape = guide_legend(override.aes = list(size = 1), order=2, reverse=TRUE),
           color = guide_legend(override.aes = list(size = 0.5), order=1)) +
    theme(legend.key.size = unit(0.5, "lines"))
  p2_species
  savePNGPDF(paste0(OUTDIR, "2-", ihh, "-", ispecies_id),
             p2_species + guides(shape="none", color="none"), height, width)
  savePNGPDF(paste0(OUTDIR, "2-speciesTrajectory-legend"),
             get_legend(p2_species), 1.25, 0.5)
}

# Plot examples of strain dynamics
plotSpeciesTrajectory("XD", "Bacteroides_xylanisolvens_57185", 1,2.1)
plotSpeciesTrajectory("XD", "Faecalibacterium_prausnitzii_62201", 1,2.1)
plotSpeciesTrajectory("XD", "Eubacterium_eligens_61678", 1,2.1)

# Plot examples of colonization and transmission.
plotSpeciesTrajectory("XD", "Clostridium_hathewayi_55515", 1, 2.1)
plotSpeciesTrajectory("XD", "Bacteroides_clarus_62282", 1,2.1)
plotSpeciesTrajectory("XD", "Alistipes_putredinis_61533", 1,2.1)
plotSpeciesTrajectory("XD", "Akkermansia_muciniphila_55290", 1,2.1)


# Export the trajectories of species of interest.
plotSpeciesTrajectoryAbxOnly <- function(ihh, ispecies_id, height, width){
  # Extract data for species of interest.
  dataSpecies <- dataSpeciesTrajectories %>%
    filter(hh==ihh, species_id==ispecies_id) %>%
    mutate(relative_abundance=ifelse(relative_abundance<1e-5, 1e-5, relative_abundance)) %>%
    mutate(timePeriod=ifelse(sample %in% samplesXmain, "main", "followup")) %>%
    mutate(subjectResponse=gsub(" ","\n", annotateSubjectResponse(subject))) %>%
    filter(subject %in% subjectsAbx,
           sample %in% samplesXmain)
  # Set a color for the abx subject based on the color of subject responses.
  abxSubject <- unique((dataSpecies %>% filter(subject %in% subjectsAbx))$subject)
  strainSharingShapes <- c(1,16)
  names(strainSharingShapes) <- c(FALSE,TRUE)
  abxStartIndex <- which(sort(unique(dataSpecies$timepoint))==29)
  abxEndIndex <- which(sort(unique(dataSpecies$timepoint))==34)
  followupStartIndex <- min(which(sort(unique(dataSpecies$timepoint))>75))
  followupEndIndex <- max(which(sort(unique(dataSpecies$timepoint))>75))
  p2_species <- dataSpecies  %>%
    arrange((subjectResponse)) %>%
    ggplot() +
    scale_y_continuous(trans='log10', limits=c(1e-5,1),
                       breaks=trans_breaks('log10', function(x) 10^x),
                       labels=trans_format('log10', math_format(10^.x))) +
    geom_hline(yintercept=1e-5, linetype="dashed") +
    annotate(geom="rect", xmin=0, xmax=5, ymin=10^-5, ymax=1,
             fill=abxColor, alpha=0.5) +
    geom_line(aes(x=(timepoint-29), y=relative_abundance, group=subject,
                  color=factor(subjectResponse))) +
    geom_point(aes(x=(timepoint-29), y=relative_abundance,
                   color=factor(subjectResponse),
                   shape=factor(relative_abundance==1e-5)), alpha=0.8, size=1) +
    scale_color_manual(values=PALETTE.ABXRESPONSELINEBREAK, name="Subject") +
    scale_shape_manual(values=c(16,1)) +
    xlab("Study day") + ylab("Relative\nabundance") +
    DEFAULTS.THEME_PRINT +
    ggtitle(gsub("_"," ",substr(ispecies_id, 1, nchar(ispecies_id)-6))) +
    theme(plot.title=element_text(size=5, margin=margin(0,0,1,0)),
          strip.text.x=element_blank()) +
    theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=5),
          axis.text.y=element_text(size=5),
          axis.title=element_text(size=7)) +
    guides(shape = guide_legend(override.aes = list(size = 1), order=2, reverse=TRUE),
           color = guide_legend(override.aes = list(size = 0.5), order=1)) +
    theme(legend.key.size = unit(0.5, "lines"))
  p2_species
  savePNGPDF(paste0(OUTDIR, "2-", ihh, "-", ispecies_id, "abxOnly"),
             p2_species + guides(shape="none", color="none"), height, width)
}

# Export the trajectories of species of interest.
plotSpeciesTrajectoryAbxOnlyFollowup <- function(ihh, ispecies_id, height, width){
  # Extract data for species of interest.
  dataSpecies <- dataSpeciesTrajectories %>%
    filter(hh==ihh, species_id==ispecies_id) %>%
    mutate(relative_abundance=ifelse(relative_abundance<1e-5, 1e-5, relative_abundance)) %>%
    mutate(timePeriod=ifelse(sample %in% samplesXmain, "main", "followup")) %>%
    mutate(subjectResponse=gsub(" ","\n", annotateSubjectResponse(subject))) %>%
    filter(subject %in% subjectsAbx)
  # Set a color for the abx subject based on the color of subject responses.
  abxSubject <- unique((dataSpecies %>% filter(subject %in% subjectsAbx))$subject)
  strainSharingShapes <- c(1,16)
  names(strainSharingShapes) <- c(FALSE,TRUE)
  abxStartIndex <- which(sort(unique(dataSpecies$timepoint))==29)
  abxEndIndex <- which(sort(unique(dataSpecies$timepoint))==34)
  followupStartIndex <- min(which(sort(unique(dataSpecies$timepoint))>75))
  followupEndIndex <- max(which(sort(unique(dataSpecies$timepoint))>75))
  p2_species <- dataSpecies  %>%
    arrange((subjectResponse)) %>%
    ggplot() +
    scale_y_continuous(trans='log10', limits=c(1e-5,1),
                       breaks=trans_breaks('log10', function(x) 10^x),
                       labels=trans_format('log10', math_format(10^.x))) +
    geom_hline(yintercept=1e-5, linetype="dashed") +
    annotate(geom="rect", xmin=0, xmax=5, ymin=10^-5, ymax=1,
             fill=abxColor, alpha=0.5) +
    geom_line(aes(x=(timepoint-29), y=relative_abundance, group=subject,
                  color=factor(subjectResponse))) +
    geom_point(aes(x=(timepoint-29), y=relative_abundance,
                   color=factor(subjectResponse),
                   shape=factor(relative_abundance==1e-5)), alpha=0.8, size=1) +
    scale_color_manual(values=PALETTE.ABXRESPONSELINEBREAK, name="Subject") +
    scale_shape_manual(values=c(16,1)) +
    xlab("Study day") + ylab("Relative\nabundance") +
    DEFAULTS.THEME_PRINT +
    ggtitle(gsub("_"," ",substr(ispecies_id, 1, nchar(ispecies_id)-6))) +
    theme(plot.title=element_text(size=5, margin=margin(0,0,1,0)),
          strip.text.x=element_blank()) +
    theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=5),
          axis.text.y=element_text(size=5),
          axis.title=element_text(size=7)) +
    guides(shape = guide_legend(override.aes = list(size = 1), order=2, reverse=TRUE),
           color = guide_legend(override.aes = list(size = 0.5), order=1)) +
    theme(legend.key.size = unit(0.5, "lines"))
  p2_species
  savePNGPDF(paste0(OUTDIR, "2-", ihh, "-", ispecies_id, "abxOnly"),
             p2_species + guides(shape="none", color="none"), height, width)
}

# Plot examples of strain dynamics
plotSpeciesTrajectoryAbxOnly("XD", "Bacteroides_xylanisolvens_57185", 1,2.1)
plotSpeciesTrajectoryAbxOnly("XD", "Faecalibacterium_prausnitzii_62201", 1,2.1)
plotSpeciesTrajectoryAbxOnly("XD", "Eubacterium_eligens_61678", 1,2.1)
plotSpeciesTrajectoryAbxOnlyFollowup("XD", "Bacteroides_clarus_62282", 1,2.1)
plotSpeciesTrajectoryAbxOnlyFollowup("XD", "Alistipes_putredinis_61533", 1,2.1)
plotSpeciesTrajectoryAbxOnlyFollowup("XD", "Parabacteroides_johnsonii_55217", 1, 2.1)
plotSpeciesTrajectoryAbxOnlyFollowup("XA", "Alistipes_putredinis_61533", 1, 2.1)


# Figure 3 - Limited colonization in the weeks after antibiotic pe --------

fig3summarywidth <- 1.25
fig3height <- 1.25
fig3relativeabundancewidth <- 1.8
fig3trajectorywidth <- 3.75


source("workflow/analysis/plotHelpers.R")
OUTDIR <- "workflow/analysis/figures-v6/out/"

# Import the JSD values for each subject from the initial timepoint.
dataBeta <- read.table("workflow/analysis/calculateDiversityStats/out/speciesBeta.txt.gz",
                         header=TRUE, stringsAsFactors = FALSE)
# Extract only JSD between the initial timepoint and
# other timepoints from the same subject
dataBetaJSDInitialSameSubject <- dataBeta %>%
  filter(method=="jsd", sample1 %in% samplesInitial,
         substr(sample1,1,3)==substr(sample2,1,3)) %>%
  mutate(subject=substr(sample1,1,3), timepoint=as.numeric(substr(sample2,5,7)))
# Calculate the maximum and final JSD From the initial timepoint.
dataBetaJSDMaxFinal <- dataBetaJSDInitialSameSubject %>%
  filter(sample2 %in% samplesXmain) %>%
  group_by(subject) %>%
  summarize(maxJSD=max(value),
            finalJSD=value[timepoint==max(timepoint)]) %>%
  mutate(subjectResponse=annotateSubjectResponse(subject))
abxSubjectsOrderedJSD <- dataBetaJSDMaxFinal %>%
  filter(subjectResponse!="control") %>%
  arrange(maxJSD) %>% pull(subject)

# Plot the number of species/strains that are maintained, disrupted, recovered, and colonized.
# Import the number of species in each category.
dataNumSpeciesPerCategory <-
  read.table("workflow/analysis/annotateSpeciesTrajectoriesSharedStrains/out/speciesTrajectoriesPerSubject-summary.txt",
             header=TRUE, stringsAsFactors = FALSE)
# Create a smaller number of categories.
dataNumSpeciesPerCategory <- dataNumSpeciesPerCategory %>%
  mutate(speciesTrajectory=ifelse(speciesTrajectory %in% c("speciesColonizedPreAbx", "speciesColonizedAbx",
                                                           "speciesColonizedPostAbx", "strainTurnover"),
                                  "colonized", speciesTrajectory)) %>%
  group_by(subject, speciesTrajectory) %>% summarize(numSpecies=sum(numSpecies)) %>%
  ungroup() %>% filter(speciesTrajectory %in% c("colonized","disrupted","recoveredMain","totalResident")) %>%
  pivot_wider(names_from=speciesTrajectory, values_from=numSpecies) %>%
  mutate(maintained=totalResident-disrupted, disruptedNotRecoveredMain=disrupted-recoveredMain) %>%
  pivot_longer(-subject, names_to="speciesTrajectory", values_to="numSpecies") %>%
  filter(speciesTrajectory %in% c("colonized","disruptedNotRecoveredMain", "recoveredMain","maintained")) %>%
  mutate(speciesTrajectory = case_when(
    speciesTrajectory=="colonized" ~ "colonized",
    speciesTrajectory=="disruptedNotRecoveredMain" ~ "notRecovered",
    speciesTrajectory=="recoveredMain" ~ "recovered",
    speciesTrajectory=="maintained" ~ "notDisrupted"
  )) %>% mutate(speciesTrajectory=fct_relevel(speciesTrajectory, c("colonized","notRecovered","recovered","notDisrupted")))
PALETTE.SPECIESTRAJECTORIES <- c("gray30","#117733","#332288","#882255")
names(PALETTE.SPECIESTRAJECTORIES) <- c("notDisrupted","recovered","notRecovered","colonized")
p3_numSpeciesPerCategory <- dataNumSpeciesPerCategory %>%
  mutate(subjectResponse=annotateSubjectResponse(subject),
         subjectResponseDisplay=annotateSubjectResponseNumSubjects(subject)) %>%
  filter(subjectResponse!="control") %>%
  ggplot() +
  geom_bar(aes(x=fct_relevel(subject, abxSubjectsOrderedJSD), y=numSpecies,
               fill=fct_relevel(speciesTrajectory, rev(names(PALETTE.SPECIESTRAJECTORIES)))),
           stat="identity") +
  facet_grid(.~fct_rev(gsub(" \\(", "\n\\(", subjectResponseDisplay)), scales="free_x", space="free_x") +
  scale_fill_manual(values=PALETTE.SPECIESTRAJECTORIES, name="Species trajectory",
                    labels=c("new species, detected in first month post-abx","species disrupted, not recovered in first month post-abx",
                             "species disrupted, recovered in first month post-abx", "species maintained")) +
  xlab("Subject") + ylab("Number of\nspecies") +
  guides(fill = guide_legend(override.aes = list(size = 1))) +
  theme(legend.key.size = unit(0.5, "lines")) +
  DEFAULTS.THEME_PRINT + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=5))
p3_numSpeciesPerCategory
savePNGPDF(paste0(OUTDIR,"3-numSpeciesPerSubject-perAnnotation"),
           p3_numSpeciesPerCategory, fig3height+0.25, 5.9)

# Calculate the percentage of disrupted species that recover in subjects with transient
# antibiotic responses.
# transient: 574 recovered / 821 disrupted, 70%
# lasting: 78 recovered / 309 disrupted, 25%
dataNumSpeciesPerCategory %>%
  mutate(subjectResponse=annotateSubjectResponse(subject)) %>%
  group_by(subjectResponse, speciesTrajectory) %>%
  summarize(numSpecies=sum(numSpecies)) %>%
  ungroup() %>% pivot_wider(names_from=speciesTrajectory, values_from=numSpecies) %>%
  mutate(pctRecovery=recovered/(recovered+notRecovered),
         totDisrupted=recovered+notRecovered)


# Import the total abundance of disrupted and colonizing species,
# annotated based only on the key timepoints.
# Import the total percentage of disrupted and colonized species at each point in time.
dataTotalDisruptedColonized <-
  read.table("workflow/analysis/annotateSpeciesTrajectoriesSharedStrains/out/speciesTrajectoriesPerSubject-totalAbundance-keyTimepointsFollowup.txt",
             header=TRUE, stringsAsFactors = FALSE)

# Extract the median and range of species disruption among control subjects.
# min: 0, median: 0.4%, max: 9%
dataTotalDisruptedColonized %>%
  mutate(subjectResponse=annotateSubjectResponse(subject)) %>%
  filter(subjectResponse=="control", speciesTrajectory=="disrupted") %>%
  group_by(subject) %>% filter(timepoint==min(timepoint)) %>% ungroup() %>%
  summarize(medianAbundanceDisrupted=median(totalAbundance),
            minAbundanceDisrupted=min(totalAbundance), maxAbundanceDisrupted=max(totalAbundance))

# Extract the median and range of species colonization and strain turnover among control subjects
# during the main study.
# min: 0, median: 0, max: 2.4%
dataTotalDisruptedColonized %>%
  mutate(subjectResponse=annotateSubjectResponse(subject)) %>%
  filter(subjectResponse=="control", speciesTrajectory %in% c("colonized", "strainTurnover"),
         sample %in% samplesFinal) %>%
  group_by(subject) %>%
  summarize(speciesTrajectory="colonizedStrainTurnover", totalAbundance=sum(totalAbundance)) %>%
  ungroup() %>%
  summarize(medianAbundanceColonized=median(totalAbundance),
            minAbundanceColonized=min(totalAbundance), maxAbundanceColonized=max(totalAbundance))
# Extract the number of control subjects in whom we detected no colonization during the main study.
# 20 subjects, 20/26, 77%
dataTotalDisruptedColonized %>%
  mutate(subjectResponse=annotateSubjectResponse(subject)) %>%
  filter(subjectResponse=="control", speciesTrajectory %in% c("colonized", "strainTurnover"),
         sample %in% samplesFinal) %>%
  group_by(subject) %>%
  summarize(speciesTrajectory="colonizedStrainTurnover", totalAbundance=sum(totalAbundance)) %>%
  ungroup() %>%
  filter(totalAbundance==0) %>% pull(subject) %>% length()
# Calculate the total number of control subjects
# There are 26 control subjects and 22 antibiotic-taking subjects.
dataTotalDisruptedColonized %>%
  mutate(subjectResponse=annotateSubjectResponse(subject)) %>%
  filter(subjectResponse=="control") %>% pull(subject) %>% unique() %>% length()


# Import the total abundance of species over time,
# annotated based on all timepoints.
dataTotalDisruptedColonized_allTimepoints <-
  read.table("workflow/analysis/annotateSpeciesTrajectoriesSharedStrains/out/speciesTrajectoriesPerSubject-totalAbundance-allTimepoints.txt",
             header=TRUE, stringsAsFactors = FALSE)

# Calculate the relative abundance of disrupted species at the beginning of the main study.
# minimal response - min: 19%, median:40%, max: 50%
# transient response - min: 73%, median: 75%, max: 79%
# lasting response - min: 69%, median: 89%, max: 96%
dataTotalDisruptedColonized_allTimepoints %>%
  mutate(subjectResponse=annotateSubjectResponse(subject)) %>%
  filter(sample %in% samplesInitial, speciesTrajectory=="disrupted") %>%
  group_by(subjectResponse) %>%
  summarize(minAbundance=min(totalAbundance), medianAbundance=median(totalAbundance),
            maxAbundance=max(totalAbundance), numZero=sum(totalAbundance==0), totalSubjects=n())
# Calculate the relative abundance of disrupted species at the end of the main study.
# minimal responses: min: 9%, median:34%, max: 50%
# transient responses: min: 43%, median: 70%, max: 78%
# lasting responses: min: 0.1%, median: 17%, max: 56%
dataTotalDisruptedColonized_allTimepoints %>%
  mutate(subjectResponse=annotateSubjectResponse(subject)) %>%
  filter(sample %in% samplesFinal, speciesTrajectory=="disrupted") %>% group_by(subjectResponse) %>%
  summarize(minAbundance=min(totalAbundance), medianAbundance=median(totalAbundance),
            maxAbundance=max(totalAbundance), numZero=sum(totalAbundance==0), totalSubjects=n())

# Calculate the median amount of species colonization at the end of the main study.
# minimal response: min: 0, median: 0, max: 9%, 5/9 subjects with no colonization detected (56%)
# transient response: min: 0, median: 0, max: 7%, 4/8 subjects with no colonization detected (50%)
# lasting response: min: 0, median: 4.9%, max: 62%, 0/5 subjects with no colonization detected (0%)
dataTotalDisruptedColonized_allTimepoints %>%
  mutate(subjectResponse=annotateSubjectResponse(subject)) %>%
  filter(sample %in% samplesFinal, speciesTrajectory %in% c("colonized", "strainTurnover")) %>%
  group_by(subject, subjectResponse) %>%
  summarize(speciesTrajectory="colonizedStrainTurnover", totalAbundance=sum(totalAbundance)) %>%
  ungroup() %>% group_by(subjectResponse) %>%
  summarize(minAbundance=min(totalAbundance), medianAbundance=median(totalAbundance),
                          maxAbundance=max(totalAbundance), numZero=sum(totalAbundance==0),
                          totalSubjects=n())

# Import the annotated species trajectories to plot the number of new colonizers over time.
dataSpeciesTrajectories <-
  read.table("workflow/analysis/annotateSpeciesTrajectoriesSharedStrains/out/speciesAbundances-trajectoriesSharingAnnotated-all.txt.gz",
             header=TRUE, stringsAsFactors = FALSE)

# Plot the total number of newly colonizing strains at each point in time.
dataColonizingSpeciesMainTimeline <- dataSpeciesTrajectories %>%
  filter(timepoint<75) %>%
  mutate(speciesTrajectory=ifelse((timeOfColonizationAllTimepoints<75 | timeOfStrainTurnover<75),
                                  "colonizedStrainTurnover", speciesTrajectory),
         speciesTrajectory=ifelse(is.na(speciesTrajectory), "none", speciesTrajectory),
         timeOfStrainColonization=ifelse(!is.na(timeOfColonizationAllTimepoints),
                                         timeOfColonizationAllTimepoints, timeOfStrainTurnover),
         timeOfStrainColonization=ifelse(is.na(timeOfStrainColonization),0,timeOfStrainColonization)) %>%
  dplyr::select(subject, timepoint, species_id, speciesTrajectory, timeOfStrainColonization) %>%
  group_by(subject, timepoint) %>%
  summarize(cumulativeColonization=sum(speciesTrajectory=="colonizedStrainTurnover" &
                                         timeOfStrainColonization<=timepoint))
# Output the total number of colonization and strain turnover events in each subject.
# XAA: 17; XBA: 11; XDA: 4; XEA: 4; XKA: 2
# total 60 colonization/strain turnover events among all subjects, 53 in antibiotic-taking subjects
dataColonizingSpeciesMainTimeline %>%
  group_by(subject) %>% filter(timepoint==max(timepoint)) %>%
  mutate(subjectResponse=annotateSubjectResponse(subject)) %>%
  filter(subjectResponse!="control") %>% group_by(subjectResponse) %>% filter(subjectResponse=="lasting response")

# Extract a list of colonization and strain turnover events during the main study.
speciesTrajectoriesSummary <-
  read.table("workflow/analysis/annotateSpeciesTrajectoriesSharedStrains/out/speciesTrajectoriesSharingSummary.txt",
             header=TRUE, stringsAsFactors = FALSE)
speciesColonizationTurnoverMain <- speciesTrajectoriesSummary %>%
  filter(subject %in% subjectsAbx,
         (speciesStatusColonizationAllTimepoints!="speciesNotColonized" & timeOfColonizationAllTimepoints<75)|
           (!is.na(strainTurnover) & strainTurnover & timeOfStrainTurnover<75)) %>%
  dplyr::select(subject, species_id, family, timeOfColonizationAllTimepoints, timeOfStrainTurnover, speciesColonizedTransientlyMain)
# Determine the proportion of species colonization events during the main study that are transient.
# transient species: 22/47 (47%); non-transient species: 25/47
speciesColonizationTurnoverMain %>%
  dplyr::count(speciesColonizedTransientlyMain)


# Plot a timeline of recovery and colonization events during the main study.
# Import the species trajectory fits.
dataTrajectoryFits <- read.table("workflow/analysis/fitColonizationTrajectories/out/speciesTrajectoriesFit-summary.txt",
                                 header=TRUE, stringsAsFactors = FALSE, sep="\t")
# Summarize the events that occur during the main study in each group of antibiotic-taking subjects.
# Exclude events that involve unknown strains from this visualization.
dataTrajectoryFitsMain <- dataTrajectoryFits %>%
  filter(tstarTimepoint<75, typeOfRecoveryColonizationTurnover!="unknown strain",
         subject %in% subjectsAbx) %>%
  mutate(subjectResponse=annotateSubjectResponse(subject)) %>%
  mutate(speciesColonizedTransientlyMain=
           ifelse(is.na(speciesColonizedTransientlyMain), FALSE, speciesColonizedTransientlyMain)) %>%
  group_by(subjectResponse, tstarTimepoint) %>%
  arrange(subjectResponse, tstarTimepoint,
          desc(typeOfRecoveryColonizationTurnover), (speciesColonizedTransientlyMain)) %>%
  mutate(timeDetectedResponseIndex=row_number())
# Summarize the events that occur during the follow-up study in each group of antibiotic-taking subjects.
# Exclude events that involve unknown strains from this visualization.
dataTrajectoryFitsFollowup <- dataTrajectoryFits %>%
  filter(tstarTimepoint>75, typeOfRecoveryColonizationTurnover!="unknown strain",
         subject %in% subjectsAbx) %>%
  mutate(subjectResponse=annotateSubjectResponse(subject)) %>%
  mutate(speciesColonizedTransientlyAll=
           ifelse(is.na(speciesColonizedTransientlyAll), FALSE, speciesColonizedTransientlyAll)) %>%
  mutate(tstarTimepoint=tstarTimepoint-29, tstarTimepoint=floor(tstarTimepoint/100)*100+100) %>%
  group_by(subjectResponse, tstarTimepoint) %>%
  arrange(subjectResponse, tstarTimepoint,
          desc(typeOfRecoveryColonizationTurnover), (speciesColonizedTransientlyAll)) %>%
  mutate(timeDetectedResponseIndex=row_number())
# Create a dataset for geom_blank with the max for each group from the main dataset.
dataTrajectoryFitsMainLimits <- dataTrajectoryFitsMain %>%
  mutate(subjectResponseDisplay=annotateSubjectResponseNumSubjects(subject)) %>%
  group_by(subjectResponseDisplay, tstarTimepoint) %>% summarize(numSpecies=n()) %>%
  ungroup() %>% group_by(subjectResponseDisplay) %>%
  summarize(maxSpecies=max(numSpecies), tstarMax=tstarTimepoint[numSpecies==max(numSpecies)])
# Create a dataset for geom_blank with the max for each group from the followup dataset.
dataTrajectoryFitsFollowupLimits <- dataTrajectoryFitsFollowup %>%
  mutate(subjectResponseDisplay=annotateSubjectResponseNumSubjects(subject)) %>%
  group_by(subjectResponseDisplay, tstarTimepoint) %>% summarize(numSpecies=n()) %>%
  ungroup() %>% group_by(subjectResponseDisplay) %>%
  summarize(maxSpecies=max(numSpecies), tstarMax=tstarTimepoint[numSpecies==max(numSpecies)])

# Plot a barplot of species recovery and colonization events in each subject response class.
p3_timelineColonizationBar <- dataTrajectoryFitsMain %>%
  mutate(typeOfRecoveryColonizationTurnover=
           ifelse(speciesColonizedTransientlyMain, "new species, transient",
           ifelse(typeOfRecoveryColonizationTurnover=="new species",
                  "new species, persistent", typeOfRecoveryColonizationTurnover))) %>%
  mutate(tstarTimepoint=tstarTimepoint-29) %>%
  mutate(subjectResponseDisplay=annotateSubjectResponseNumSubjects(subject)) %>%
  ggplot() +
  geom_blank(data=dataTrajectoryFitsFollowupLimits, aes(x=(1), y=maxSpecies)) +
  geom_rect(xmin=0, xmax=5, ymin=0, ymax=50, fill=abxColor) +
  geom_bar(aes(x=(tstarTimepoint),
               fill=fct_relevel(typeOfRecoveryColonizationTurnover, rev(names(PALETTESTRAINS))),
               color=fct_relevel(typeOfRecoveryColonizationTurnover, rev(names(PALETTESTRAINS))))) +
  facet_wrap(~fct_rev(subjectResponseDisplay), scales="free", ncol=1) +
  scale_fill_manual(values=PALETTESTRAINSFILL, name="") +
  scale_color_manual(values=PALETTESTRAINS, name="") +
  scale_alpha_identity() + scale_y_continuous(breaks=pretty_breaks()) +
  theme(legend.key.size = unit(0.5, "lines")) +
  xlab("Colonization/recovery time") + ylab("Number of events") +
  guides(fill = guide_legend(override.aes = list(size = 0.5)), color="none") +
  DEFAULTS.THEME_PRINT + theme(legend.position="bottom")
p3_timelineColonizationBar
savePNGPDF(paste0(OUTDIR, "3-timeline-recoveryColonization-bar"),
           p3_timelineColonizationBar+ guides(fill="none"), 2.5, 2.3)
savePNGPDF(paste0(OUTDIR, "3-timeline-recoveryColonization-bar-legend"),
           get_legend(p3_timelineColonizationBar), 0.5, 3.5)


# Plot a barplot of species recovery and colonization events in each subject response class.
p3_timelineColonizationFollowupBar <- dataTrajectoryFitsFollowup %>%
  mutate(typeOfRecoveryColonizationTurnover=
           ifelse(speciesColonizedTransientlyAll, "new species, transient",
                  ifelse(typeOfRecoveryColonizationTurnover=="new species",
                         "new species, persistent", typeOfRecoveryColonizationTurnover))) %>%
  mutate(subjectResponseDisplay=annotateSubjectResponseNumSubjects(subject)) %>%
  ggplot() +
  geom_blank(data=dataTrajectoryFitsMainLimits, aes(x=(tstarMax), y=maxSpecies)) +
  geom_bar(aes(x=(tstarTimepoint),
               fill=fct_relevel(typeOfRecoveryColonizationTurnover, rev(names(PALETTESTRAINS))),
               color=fct_relevel(typeOfRecoveryColonizationTurnover, rev(names(PALETTESTRAINS))))) +
  facet_wrap(~fct_rev(subjectResponseDisplay),  ncol=1, scales="free") +
  scale_fill_manual(values=PALETTESTRAINSFILL, name="") +
  scale_color_manual(values=PALETTESTRAINS, name="") +
  scale_alpha_identity() + scale_y_continuous(breaks=pretty_breaks()) +
  scale_x_continuous(limits=c(40,1000), breaks=seq(250,1000,250)) +
  theme(legend.key.size = unit(0.5, "lines")) +
  xlab("Time detected (study day)") + ylab("") +
  guides(fill = guide_legend(override.aes = list(size = 0.5)), color="none") +
  DEFAULTS.THEME_PRINT + theme(legend.position="bottom") + theme(axis.title.y=element_blank()) +
  DEFAULTS.THEME_NOYAXIS
p3_timelineColonizationFollowupBar
savePNGPDF(paste0(OUTDIR, "3-timeline-recoveryColonization-followup-bar"),
           p3_timelineColonizationFollowupBar + guides(fill="none"), 2.5, 1.25)

# Calculate the proportion of species recoveries in subjects with minimal responses
# that are driven by pre-existing strains versus the colonization of new strains.
# minimal response: 6 new species, 1 new strain, 42 same strain, 68 unknown strain (111 strain events)
# transient response: 5 new species, 3 new strain, 87 same strain, 85 unknown strain (175 strain events)
# lasting response: 35 new species, 4 new strain, 14 same strain, 37 unknown strain (55 strain events)
# minimal+transient: 129/286, 45% resident strain; 4/286, 1.3% new strain; 153/286, 53% strain unknown
dataTrajectoryFits %>%
  mutate(subjectResponse=annotateSubjectResponse(subject)) %>%
  filter(tstarTimepoint<75) %>%
  group_by(subjectResponse, typeOfRecoveryColonizationTurnover) %>% summarize(numEvents=n())

# Identify the cases in which resident strains do not recover after antibiotics
# and we detect a new strain instead. This is only subjects with transient responses.
dataNewStrainsTransientRespondersMain <- dataTrajectoryFits %>%
  filter(annotateSubjectResponse(subject)=="transient response",
         typeOfRecoveryColonizationTurnover=="new strain",
         timeOfRecoveryColonizationTurnover<75) %>%
  dplyr::select(subject, species_id, timeOfRecoveryColonizationTurnover)
# Import the species trajectory annotations and determine if these
# new strains are shared with cohabiting partners after the colonization.
# Import the annotated species trajectories.
dataSpeciesTrajectories <- 
  read.table("workflow/analysis/annotateSpeciesTrajectoriesSharedStrains/out/speciesAbundances-trajectoriesSharingAnnotated-all.txt.gz",
             header=TRUE, stringsAsFactors = FALSE)
# Pare down the dataframe so that it only contains information about strain sharing
# across various timepoints.
dataSpeciesTrajectoriesSharing <- dataSpeciesTrajectories %>%
  dplyr::select(subject, timepoint, species_id, relative_abundance, sharingHhAnyTimepoint)
# Extract the abundances of new strains and determine if they were shared.
# Of the four species with new strains, three are not shared, and the strain sharing
# status of the fourth species is unknown.
left_join(dataNewStrainsTransientRespondersMain,
          dataSpeciesTrajectoriesSharing, by=c("subject","species_id")) %>%
  filter(timepoint>=timeOfRecoveryColonizationTurnover)


# Identify all cases of colonization for which we successfully fit a trajectory.
# Determine if each of these colonization events was due to transmission
# from the cohabiting partner.
dataNewColonizers <- dataTrajectoryFits %>%
  filter(typeOfRecoveryColonizationTurnover %in% c("new strain", "new species")) %>%
  dplyr::select(subject, species_id, timeOfRecoveryColonizationTurnover,
                speciesColonizedTransientlyAll)
# Extract the abundances of new colonizers and determine if they were shared
# with the cohabiting partner.
dataNewColonizersSharing <- left_join(dataNewColonizers, dataSpeciesTrajectoriesSharing, 
          by=c("subject","species_id")) %>%
  filter(timepoint>=timeOfRecoveryColonizationTurnover)
# Summarize the strain sharing calls for each species.
dataNewColonizersSharingSummary <- dataNewColonizersSharing %>%
  group_by(subject, species_id, speciesColonizedTransientlyAll, sharingHhAnyTimepoint) %>%
  summarize(numTimepoints=n(), finalAbundance=relative_abundance[timepoint==max(timepoint)]) %>%
  mutate(transient=(finalAbundance<1e-3)) %>%
  pivot_wider(names_from=sharingHhAnyTimepoint, values_from=numTimepoints,
              values_fill=0) %>%
  mutate(totalTimepoints=strainsUnknown+strainsShared+speciesNotShared+strainsNotShared) %>%
  mutate(sharingStatus=case_when(
    strainsUnknown==totalTimepoints ~ "strainsUnknown",
    strainsShared==totalTimepoints ~ "strainsShared",
    speciesNotShared==totalTimepoints ~ "speciesNotShared",
    strainsNotShared==totalTimepoints ~ "strainsNotShared"
  ))
# Summarize the origin of new colonizers, across all colonizers.
# Total: 149 colonization events (79 with strain status known)
# species not shared: 43/149 (29%) of all colonization events, 43/79 (54%) of colonization events with known status
# strains not shared: 24/149 (16%) of all colonization events, 24/79 (30%) of colonization events with known status
# strains shared: 12/79 (8%) of all colonization events, 12/79 (15%) of all colonization events
dataNewColonizersSharingSummary %>% ungroup() %>%
  count(sharingStatus) %>%
  mutate(pctStatus=n/sum(n))
# Summarize the origin of new colonizers for colonizers that are still detected
# at the end of the study in all subjects.
# Total: 81 colonization events (52 with strain status known, 29 (36%) with strain status unknown)
# species not shared: 19/81 (24%) of all colonization events, 19/52 (37%) of colonization events with known status
# strains not shared: 21/81 (26%) of all colonization events, 21/52 (40%) of colonization events with known status
# strains shared: 12/81 (15%) of all colonization events, 12/52 (23%) of all colonization events
dataNewColonizersSharingSummary %>% ungroup() %>%
  filter(!transient) %>%
  count(sharingStatus) %>%
  mutate(pctStatus=n/sum(n))

# Summarize the origin of new colonizers for colonizers that are still detected
# at the end of the study only in subjects with lasting responses.
# Total: 36 colonization events (25 with strain status known, 11 (36%) with strain status unknown)
# species not shared: 5/36 (31%) of all colonization events, 5/25 (20%) of colonization events with known status
# strains not shared: 12/36 (33%) of all colonization events, 12/25 (48%) of colonization events with known status
# strains shared: 8/36 (22%) of all colonization events, 8/25 (32%) of all colonization events
dataNewColonizersSharingSummary %>% ungroup() %>%
  filter(!transient, annotateSubjectResponse(subject)=="lasting response") %>%
  count(sharingStatus) %>%
  mutate(pctStatus=n/sum(n))


# Import the HMP prevalence and abundance of resident species and new colonizers.
residentSpeciesColonizersHMPprevalenceAbundance <-
  read.table("workflow/analysis/calculateHMPspeciesPrevalenceAbundance/out/residentColonizingSpecies-prevalenceAbundance.txt",
             header=TRUE, stringsAsFactors = FALSE) %>%
  mutate(type=fct_relevel(type, c("resident\nspecies", "new colonizer\n<1 month post-abx", "new colonizer\n>1 month post-abx")))
p3_prevalenceResidentColonizers <- residentSpeciesColonizersHMPprevalenceAbundance %>%
  ggplot() +
  geom_violin(aes(x=type, y=prevalence), fill="gray80") +
  geom_boxplot(aes(x=type, y=prevalence), width=0.1) +
  ylab("Species\nprevalence, HMP") +
  DEFAULTS.THEME_PRINT + theme(axis.title.x=element_blank())
# savePNGPDF(paste0(OUTDIR, "3-HMPprevalenceAbundance-residentsColonizers"),
#            p3_prevalenceResidentColonizers + p3_medianAbundanceResidentColonizers, fig3height, 3)
savePNGPDF(paste0(OUTDIR, "3-HMPprevalenceAbundance-residentsColonizers"),
           p3_prevalenceResidentColonizers, fig3height-0.25, 3.2)


# Compare the distributions of species prevalence between categories using the Wilcoxon rank-sum test.
# resident vs initial month colonizers: 5.01e-16
# resident vs follow-up colonizers: 2.702e-9
# initial month vs follow-up colonizers: 0.001
wilcox.test(residentSpeciesColonizersHMPprevalenceAbundance %>% filter(type=="resident\nspecies") %>% pull(prevalence),
            residentSpeciesColonizersHMPprevalenceAbundance %>% filter(type=="new colonizer\n<1 month post-abx") %>% pull(prevalence))
wilcox.test(residentSpeciesColonizersHMPprevalenceAbundance %>% filter(type=="resident\nspecies") %>% pull(prevalence),
            residentSpeciesColonizersHMPprevalenceAbundance %>% filter(type=="new colonizer\n>1 month post-abx") %>% pull(prevalence))
wilcox.test(residentSpeciesColonizersHMPprevalenceAbundance %>% filter(type=="new colonizer\n<1 month post-abx") %>% pull(prevalence),
            residentSpeciesColonizersHMPprevalenceAbundance %>% filter(type=="new colonizer\n>1 month post-abx") %>% pull(prevalence))
# Compare the distributions of median species abundance between categories using the Wilcoxon rank-sum test.
# Use log10 species abundances.
# resident vs initial month colonizers: 3.12e-12
# resident vs follow-up colonizers: 0.0002321
# initial month vs follow-up colonizers: 0.0003161
wilcox.test(log10(residentSpeciesColonizersHMPprevalenceAbundance %>% filter(type=="resident\nspecies") %>% pull(medianAbundance)),
            log10(residentSpeciesColonizersHMPprevalenceAbundance %>% filter(type=="new colonizer\n<1 month post-abx") %>% pull(medianAbundance)))
wilcox.test(log10(residentSpeciesColonizersHMPprevalenceAbundance %>% filter(type=="resident\nspecies") %>% pull(medianAbundance)),
            log10(residentSpeciesColonizersHMPprevalenceAbundance %>% filter(type=="new colonizer\n>1 month post-abx") %>% pull(medianAbundance)))
wilcox.test(log10(residentSpeciesColonizersHMPprevalenceAbundance %>% filter(type=="new colonizer\n<1 month post-abx") %>% pull(medianAbundance)),
            log10(residentSpeciesColonizersHMPprevalenceAbundance %>% filter(type=="new colonizer\n>1 month post-abx") %>% pull(medianAbundance)))


# Import species trajectory summary data.
dataSpeciesTrajectoriesSummary <-
  read.table("workflow/analysis/annotateSpeciesTrajectoriesSharedStrains/out/speciesTrajectoriesSharingSummary.txt",
             header=TRUE, stringsAsFactors = FALSE)
# Plot the relationship of strain/species sharing between subjects
# and the dynamics of recovery.
dataSharingRecovery <- dataSpeciesTrajectoriesSummary %>%
  filter(subject %in% subjectsAbx,
         speciesStatusAbxAllTimepoints %in% c("speciesNotRecovered", "speciesRecoveredMain", "speciesRecoveredFollowup")) %>%
  dplyr::select(subject, speciesStatusAbxAllTimepoints, sharingHhPreAbx)
dataSharingRecovery <- dataSharingRecovery %>%
  mutate(speciesSharing=ifelse(sharingHhPreAbx=="speciesNotShared",FALSE,TRUE),
         strainSharing=case_when(
           sharingHhPreAbx=="strainsNotShared" ~ FALSE,
           sharingHhPreAbx=="strainsShared" ~ TRUE,
           TRUE ~ NA),
         speciesRecoveredMain=
           ifelse(speciesStatusAbxAllTimepoints=="speciesRecoveredMain",
                  "recovered","not recovered"))
# Summarize the relationship between strain sharing and recovery during the main study.
dataSharingRecoveryStrains <- dataSharingRecovery %>%
  group_by(speciesRecoveredMain, strainSharing) %>%
  summarize(numSpecies=n()) %>%
  filter(!is.na(strainSharing)) %>%
  group_by(strainSharing) %>%
  mutate(pctSpecies=numSpecies/sum(numSpecies),
         totalSpecies=sum(numSpecies))
sharingRecoveryStrains_chisq <-
  chisq.test(matrix(dataSharingRecoveryStrains$numSpecies, nrow=2, ncol=2))

# Calculate the relationship between strain sharing and species recovery
# for only the species in subjects with minimal and transient responses.
# chi-square p=0.57
dataSharingRecoveryMinimalTransient <- dataSpeciesTrajectoriesSummary %>%
  mutate(subjectResponse=annotateSubjectResponse(subject)) %>%
  filter(subjectResponse %in% c("minimal response", "transient response"),
         speciesStatusAbxAllTimepoints %in% c("speciesNotRecovered", "speciesRecoveredMain", "speciesRecoveredFollowup")) %>%
  dplyr::select(subject, speciesStatusAbxAllTimepoints, sharingHhPreAbx)
dataSharingRecoveryMinimalTransient <- dataSharingRecoveryMinimalTransient %>%
  mutate(speciesSharing=ifelse(sharingHhPreAbx=="speciesNotShared",FALSE,TRUE),
         strainSharing=case_when(
           sharingHhPreAbx=="strainsNotShared" ~ FALSE,
           sharingHhPreAbx=="strainsShared" ~ TRUE,
           TRUE ~ NA),
         speciesRecoveredMain=
           ifelse(speciesStatusAbxAllTimepoints=="speciesRecoveredMain",
                  "recovered","not recovered"))
# Summarize the relationship between strain sharing and recovery during the main study.
dataSharingRecoveryStrainsMinimalTransient <- dataSharingRecoveryMinimalTransient %>%
  group_by(speciesRecoveredMain, strainSharing) %>%
  summarize(numSpecies=n()) %>%
  filter(!is.na(strainSharing)) %>%
  group_by(strainSharing) %>%
  mutate(pctSpecies=numSpecies/sum(numSpecies),
         totalSpecies=sum(numSpecies))
sharingRecoveryStrainsMinimalTransient_chisq <-
  chisq.test(matrix(dataSharingRecoveryStrainsMinimalTransient$numSpecies, nrow=2, ncol=2))

# Import the data on taxonomic enrichment among disrupted species.
dataEnrichmentDisrupted <-
  read.table("workflow/analysis/calculateTaxonomicEnrichment/out/disruptedSpecies-randomDistributionSummary.txt",
             header=TRUE, stringsAsFactors = FALSE, sep="\t")
# Plot a volcano plot of the family enrichment analysis.
dataEnrichmentDisrupted <- dataEnrichmentDisrupted %>%
  mutate(log2Enrichment=log2((meanActual+1)/(meanRandom+1)),
         significant=(pvalue<lowerThreshold),
         pvalue=ifelse(pvalue==0, 1/numCounts, pvalue),
         Bonferronipvalue=pvalue*n())
p3_enrichment_disrupted <- dataEnrichmentDisrupted %>%
  ggplot() +
  geom_hline(aes(yintercept=-log10(lowerThreshold)), linetype="dashed", alpha=0.5) +
  geom_point(aes(x=log2Enrichment, y=-log10(pvalue),
                 color=factor(ifelse(!significant, "ns",
                                     ifelse(log2Enrichment>0,"sig-plus","sig-minus")))),
             alpha=0.8, size=1) +
  scale_color_manual(values=c("gray40","dodgerblue3","firebrick")) +
  guides(color="none") + ylim(0,4.5) + xlim(-1,1) +
  xlab("Fold enrichment\namong disrupted species") +
  ylab("Significance\n-log10(p-value)") +
  DEFAULTS.THEME_PRINT
savePNGPDF(paste0(OUTDIR, "3-enrichment-volcano-disrupted"), p3_enrichment_disrupted,
           fig3height, 7-fig3summarywidth-fig3trajectorywidth-0.5)

# Import the data on taxonomic enrichment among colonizing species.
dataEnrichmentColonizingTurnover <-
  read.table("workflow/analysis/calculateTaxonomicEnrichment/out/colonizingSpeciesStrainTurnover-main-randomDistributionSummary.txt",
             header=TRUE, stringsAsFactors = FALSE, sep="\t")
# Plot a volcano plot of the family enrichment analysis.
dataEnrichmentColonizingTurnover <- dataEnrichmentColonizingTurnover %>%
  mutate(log2Enrichment=log2((meanActual+1)/(meanRandom+1)),
         significant=(pvalue<lowerThreshold),
         pvalue=ifelse(pvalue==0, 1/numCounts, pvalue),
         Bonferronipvalue=pvalue*n())
p3_enrichment_colonizingTurnover <- dataEnrichmentColonizingTurnover %>%
  ggplot() +
  geom_hline(aes(yintercept=-log10(lowerThreshold)), linetype="dashed", alpha=0.5) +
  geom_point(aes(x=log2Enrichment, y=-log10(pvalue),
                 color=factor(significant)), alpha=0.8, size=1) +
  scale_color_manual(values=c("gray40","firebrick")) +
  guides(color="none") + ylim(0,4.5) + xlim(-2.5,2.5) +
  xlab("Fold enrichment\namong new colonizers") +
  ylab("Significance\n-log10(p-value)") +
  DEFAULTS.THEME_PRINT
savePNGPDF(paste0(OUTDIR, "3-enrichment-volcano-colonizingTurnover-main"), p3_enrichment_colonizingTurnover,
           fig3height, 7-fig3summarywidth-fig3trajectorywidth-0.5)


# Figure 4 - Delayed colonization after antibiotic perturbation -----------

fig4summarywidth <- 1.25
fig4height <- 1.25
fig4relativeabundancewidth <- 1.8
fig4trajectorywidth <- 3.75

source("workflow/analysis/plotHelpers.R")
OUTDIR <- "workflow/analysis/figures-v6/out/"

# Import the total abundance of species over time,
# annotated based on all timepoints.
dataTotalDisruptedColonized_allTimepoints <-
  read.table("workflow/analysis/annotateSpeciesTrajectoriesSharedStrains/out/speciesTrajectoriesPerSubject-totalAbundance-allTimepoints.txt",
             header=TRUE, stringsAsFactors = FALSE)

# Plot the total amount of colonization and strain turnover at the final sampling timepoint.
dataTotalColonizedSummary <- dataTotalDisruptedColonized_allTimepoints %>%
  filter(sample %in% samplesX, speciesTrajectory %in% c("colonized", "strainTurnover")) %>%
  group_by(subject, timepoint) %>%
  summarize(speciesTrajectory="colonizedStrainTurnover", totalAbundance=sum(totalAbundance)) %>%
  mutate(subjectResponse=annotateSubjectResponse(subject),
         subjectTreatment=annotateSubjectTreatment(subject)) %>%
  mutate(totalAbundance=ifelse(totalAbundance<1e-5, 1e-5, totalAbundance)) %>%
  mutate(hh=substr(subject,1,2), sample=paste0(subject,"-",formatC(timepoint, width=3, flag="0")))
dataTotalColonizedSummaryEndFollowup <- dataTotalColonizedSummary %>%
  ungroup() %>% group_by(subject) %>%
  filter(timepoint==max(timepoint))
p4_pctColonized_followup <- dataTotalColonizedSummary %>% ungroup() %>% group_by(subject) %>%
  filter(timepoint==max(timepoint)) %>% arrange(desc(subjectResponse)) %>%
  ggplot() +
  geom_point(aes(x=fct_rev(subjectTreatment), y=totalAbundance, color=factor(subjectResponse)),
             position=position_jitter(height=0, width=0.2, seed=0), alpha=0.8, size=1) +
  scale_color_manual(values=PALETTE.ABXRESPONSE) +
  guides(color="none") + ylim(0,1) +
  xlab("Subject") + ylab("Total abundance,\nnew colonizers") +
  DEFAULTS.THEME_PRINT
savePNGPDF(paste0(OUTDIR,"4-pctColonized-followup"), p4_pctColonized_followup,
           fig4height, fig4summarywidth)
# Calculate the range of colonization at the end of follow-up sampling.
# control: min: 0, median: 0.2%, max: 11%
# minimal + transient response: min: 0, median: 0.5%, max: 17.5%
# lasting response: min: 3%, median: 22.9%, max: 67.5%
# Wilcoxon, control vs minimal/transient: p=0.36 (problem with ties) (indistinguishable)
# Wilcoxon, control vs lasting: p=0.00095 (problem with ties) (distinguishable)
# Wilcoxon, transient vs lasting: p=0.003649 (problem with ties) (distinguishable)
dataTotalColonizedSummaryEndFollowup %>%
  group_by(subjectResponse) %>%
  summarize(minTotalAbundance=min(totalAbundance), medianTotalAbundance=median(totalAbundance),
            maxTotalAbundance=max(totalAbundance))
wilcox.test(dataTotalColonizedSummaryEndFollowup %>% filter(subjectResponse=="control") %>% pull(totalAbundance),
       dataTotalColonizedSummaryEndFollowup %>% filter(subjectResponse=="transient response") %>% pull(totalAbundance))
wilcox.test(dataTotalColonizedSummaryEndFollowup %>% filter(subjectResponse=="control") %>% pull(totalAbundance),
       dataTotalColonizedSummaryEndFollowup %>% filter(subjectResponse=="lasting response") %>% pull(totalAbundance))
wilcox.test(dataTotalColonizedSummaryEndFollowup %>% filter(subjectResponse=="transient response") %>% pull(totalAbundance),
            dataTotalColonizedSummaryEndFollowup %>% filter(subjectResponse=="lasting response") %>% pull(totalAbundance))
dataTotalColonizedSummaryEndFollowup %>% ungroup() %>%
  filter(subjectResponse %in% c("minimal response","transient response")) %>%
  summarize(minTotalAbundance=min(totalAbundance), medianTotalAbundance=median(totalAbundance),
            maxTotalAbundance=max(totalAbundance))
# Calculate the total relative abundance of colonizers at the end of follow-up sampling
# in subjects with lasting antibiotic responses and their cohabiting controls.
# XAA: 58%, XAC: 0, XAD: 0
# XBA: 3.1%, XBB, 0
# XDA: 15.2%, XDB, 3.3%
# XEA: 23%, XEB: 0.6%
# XKA: 68%, XDB: 2.3%
# abx: min: 3.1%, median: 23%, max: 68%
# control: min: 0%, median: 0.4%, max: 3.3%
dataTotalColonizedSummaryEndFollowup %>%
  filter(subject %in% c("XAA","XAC","XAD","XBA","XBB","XDA","XDB","XEA","XEB","XKA","XKB"))
dataTotalColonizedSummaryEndFollowup %>%
  filter(subject %in% c("XAA","XAC","XAD","XBA","XBB","XDA","XDB","XEA","XEB","XKA","XKB")) %>%
  group_by(subjectTreatment) %>%
  summarize(minTotalAbundance=min(totalAbundance), medianTotalAbundance=median(totalAbundance),
            maxTotalAbundance=max(totalAbundance))


# Plot the total amount of colonization and strain turnover over time.
p4_pctColonized_byResponse <- dataTotalColonizedSummary %>%
  filter(!(sample %in% samplesXBAextraFollowup)) %>%
  mutate(subjectResponseDisplay=annotateSubjectResponseNumSubjects(subject),
         timepoint=timepoint-29) %>%
  ggplot() +
  geom_line(aes(x=timepoint, y=totalAbundance, group=subject,
                color=factor(subjectResponse)),
            linewidth=0.5, alpha=0.8) +
  scale_color_manual(values=PALETTE.ABXRESPONSE, name="Subject response") +
  xlab("Study day") + ylab("Total abundance,\nnew colonizers") +
  facet_wrap(~fct_relevel(subjectResponseDisplay, subjectResponseNumSubjectsOrder), ncol=4) +
  ylim(0,1) + guides(color="none") +
  DEFAULTS.THEME_PRINT
savePNGPDF(paste0(OUTDIR,"4-pctColonizedTurnover-byResponse"), p4_pctColonized_byResponse,
           fig4height, 7-2.6)


# Import the annotated species trajectories.
dataSpeciesTrajectories <- read.table("workflow/analysis/annotateSpeciesTrajectoriesSharedStrains/out/speciesAbundances-trajectoriesSharingAnnotated-all.txt.gz",
                                      header=TRUE, stringsAsFactors = FALSE)
# Plot the cumulative number of colonizing strains during the full study.
# Include strains that colonized during the main part of the study.
dataColonizingSpeciesAllSamples <- dataSpeciesTrajectories %>%
  mutate(speciesTrajectory=ifelse((!is.na(timeOfColonizationAllTimepoints) | !is.na(timeOfStrainTurnover)),
                                  "colonizedStrainTurnover", "none"),
         speciesTrajectory=ifelse(is.na(speciesTrajectory), "none", speciesTrajectory),
         timeOfStrainColonization=ifelse(!is.na(timeOfColonizationAllTimepoints),
                                         timeOfColonizationAllTimepoints, timeOfStrainTurnover),
         timeOfStrainColonization=ifelse(is.na(timeOfStrainColonization),0,timeOfStrainColonization)) %>%
  dplyr::select(subject, timepoint, species_id, speciesTrajectory, timeOfStrainColonization) %>%
  group_by(subject, timepoint) %>%
  summarize(cumulativeColonization=sum(speciesTrajectory=="colonizedStrainTurnover" &
                                         timeOfStrainColonization<=timepoint))

# Calculate the range of the total number of colonization events identified by the end
# of follow-up sampling.
# control: min: 0, median: 1, max: 6; 9/26 (35%) no colonization
# minimal: min: 0, median: 1, max: 5; 1/9 (9%) no colonization
# transient: min: 0, median: 1.5, max: 4; 2/8 (25%) no colonization
# minimal/transient: min: 0, median: 1, max: 5; 3/17 (18%) no colonization
# lasting: min: 6, median: 11, max: 36; 0/5 (0%) no colonization
dataColonizingSpeciesAllSamples %>%
  group_by(subject) %>% filter(timepoint==max(timepoint)) %>%
  mutate(subjectResponse=annotateSubjectResponse(subject)) %>%
  ungroup() %>% group_by(subjectResponse) %>%
  summarize(minNumColonizers=min(cumulativeColonization), medianNumColonizers=median(cumulativeColonization),
            maxNumColonizers=max(cumulativeColonization),
            numZeroColonizers=sum(cumulativeColonization==0), numSubjects=n())

# Plot the number of colonizing strains during follow-up sampling.
# Do not include strains that colonized during the main part of the study.
dataColonizingSpeciesFollowup <- dataSpeciesTrajectories %>%
  mutate(speciesTrajectory=ifelse(((timeOfColonizationAllTimepoints>75) | (timeOfStrainTurnover>75)),
                                  "colonizedStrainTurnover", "none"),
         speciesTrajectory=ifelse(is.na(speciesTrajectory), "none", speciesTrajectory),
         timeOfStrainColonization=ifelse(!is.na(timeOfColonizationAllTimepoints),
                                         timeOfColonizationAllTimepoints, timeOfStrainTurnover),
         timeOfStrainColonization=ifelse(is.na(timeOfStrainColonization),0,timeOfStrainColonization)) %>%
  dplyr::select(subject, timepoint, species_id, speciesTrajectory, timeOfStrainColonization) %>%
  group_by(subject, timepoint) %>%
  summarize(cumulativeColonization=sum(speciesTrajectory=="colonizedStrainTurnover" &
                                         timeOfStrainColonization<=timepoint))

# Calculate the range of the total number of colonization events identified
# in follow-up sampling only.
# control: min: 0, median: 1, max: 6; 11/26 (42%) no colonization
# minimal: min: 0, median: 0, max: 3; 5/9 (9%) no colonization
# transient: min: 0, median: 0.5, max: 2; 4/8 (50%) no colonization
# minimal + transient: min: 0, median: 0, max: 3; 9/17 (53%) no colonization
# lasting: min: 4, median: 7, max: 19; 0/5 (0%) no colonization
# t-test: control vs minimal: p=0.26 (indistinguishable)
# t-test: control vs transient: p=0.09 (indistinguishable)
# t-test: control vs minimal+transient: p=0.11 (indistinguishable)
# t-test: control vs lasting: p=0.04 (distinguishable)
dataColonizingSpeciesFollowup %>%
  group_by(subject) %>% filter(timepoint==max(timepoint)) %>%
  mutate(subjectResponse=annotateSubjectResponse(subject)) %>%
  ungroup() %>% group_by(subjectResponse) %>%
  summarize(minNumColonizers=min(cumulativeColonization), medianNumColonizers=median(cumulativeColonization),
            maxNumColonizers=max(cumulativeColonization),
            numZeroColonizers=sum(cumulativeColonization==0), numSubjects=n())
numColonizingSpeciesPerSubjectFollowup <- dataColonizingSpeciesFollowup %>%
  group_by(subject) %>% filter(timepoint==max(timepoint)) %>%
  mutate(subjectResponse=annotateSubjectResponse(subject)) %>% ungroup()
numColonizingSpeciesPerSubjectFollowup %>%
  filter(subjectResponse %in% c("minimal response", "transient response")) %>%
  summarize(minNumColonizers=min(cumulativeColonization), medianNumColonizers=median(cumulativeColonization),
            maxNumColonizers=max(cumulativeColonization),
            numZeroColonizers=sum(cumulativeColonization==0), numSubjects=n())
t.test(numColonizingSpeciesPerSubjectFollowup %>% filter(subjectResponse=="control") %>% pull(cumulativeColonization),
       numColonizingSpeciesPerSubjectFollowup %>% filter(subjectResponse=="minimal response") %>% pull(cumulativeColonization))
t.test(numColonizingSpeciesPerSubjectFollowup %>% filter(subjectResponse=="control") %>% pull(cumulativeColonization),
       numColonizingSpeciesPerSubjectFollowup %>% filter(subjectResponse=="transient response") %>% pull(cumulativeColonization))
t.test(numColonizingSpeciesPerSubjectFollowup %>% filter(subjectResponse=="control") %>% pull(cumulativeColonization),
       numColonizingSpeciesPerSubjectFollowup %>% filter(subjectResponse %in% c("minimal response", "transient response")) %>% pull(cumulativeColonization))
t.test(numColonizingSpeciesPerSubjectFollowup %>% filter(subjectResponse=="control") %>% pull(cumulativeColonization),
       numColonizingSpeciesPerSubjectFollowup %>% filter(subjectResponse=="lasting response") %>% pull(cumulativeColonization))
# Calculate the number of follow-up sampling events in each of the subjects
# with lasting responses and their cohabiting controls.
# XAA: 19, XAC: 0, XAD: 0
# XBA: 18, XBB: 2
# XDA: 7, XDB: 2
# XEA: 7, XEB: 3
# XKA: 4, XKB: 2
dataColonizingSpeciesFollowup %>%
  group_by(subject) %>% filter(timepoint==max(timepoint)) %>%
  mutate(subjectResponse=annotateSubjectResponse(subject)) %>%
  filter(subject %in% c("XAA","XAC","XAD","XBA","XBB","XDA","XDB","XEA","XEB","XKA","XKB"))
# Calculate the number of follow-up sampling events that occur during the second half
# of follow-up sampling, where the second half is defined both by time and by the
# number of follow-up samples.
dataTimingColonizingSpeciesFollowup <- dataColonizingSpeciesFollowup %>%
  filter(timepoint>75) %>%
  mutate(subjectResponse=annotateSubjectResponse(subject)) %>%
  group_by(subject) %>% filter(subjectResponse=="lasting response") %>%
  mutate(midTimepoint=64+(max(timepoint)-64)/2, midTimepointSample=timepoint[timepoint==min(timepoint[timepoint>midTimepoint])],
         sampleIndex=row_number(), midSample=timepoint[sampleIndex==ceiling(max(sampleIndex/2))],
         previousTimepoint=lag(timepoint), previousTimepoint=ifelse(is.na(previousTimepoint),64,previousTimepoint)) %>%
  ungroup() %>% group_by(subject, midTimepointSample, midSample) %>%
  summarize(numColonizationAfterMidTimepointSample=max(cumulativeColonization)-cumulativeColonization[timepoint==midTimepointSample],
         numColonizationAfterMidSample=max(cumulativeColonization)-cumulativeColonization[timepoint==midSample],
         midColonizationTimepoint=previousTimepoint[timepoint==min(timepoint[cumulativeColonization>max(cumulativeColonization)/2])],
         totalColonizationEventsFollowup=max(cumulativeColonization))

# Import species trajectory data.
dataSpeciesTrajectoriesFit <- read.table("workflow/analysis/fitColonizationTrajectories/out/speciesTrajectoriesFit.txt",
                                         header=TRUE, stringsAsFactors = FALSE, sep="\t")
plotSpeciesTrajectory <- function(isubject, ispecies_id, height, width){
  # Export the trajectories of species of interest.
  # Extract data for species of interest.
  dataSpecies <- dataSpeciesTrajectoriesFit %>%
    filter(subject==isubject, species_id==ispecies_id) %>%
    mutate(relative_abundance=ifelse(relative_abundance<1e-5, 1e-5, relative_abundance)) %>%
    mutate(timePeriod=ifelse(sample %in% samplesXmain, "main", "followup")) %>%
    mutate(subjectResponse=gsub(" ","\n", annotateSubjectResponse(subject)))
  inferred_tau <- unique(dataSpecies$inferred_tau)
  lastUndetectableTimepoint <- unique(dataSpecies$lastUndetectableTimepointBeforeTstar)
  tstarTimepoint <- unique(dataSpecies$tstarTimepoint)
  K <- unique(dataSpecies$K)
  p4_species <- dataSpecies %>%
    mutate(timePeriod=ifelse(timepoint<75, "main", 
                             ifelse(timepoint>=lastUndetectableTimepointBeforeTstar, "new","waiting")),
           timepoint=timepoint-29) %>%
    filter(relative_abundance>1e-5 | lag(relative_abundance)>1e-5 | lead(relative_abundance)>1e-5) %>%
    ggplot() +
    geom_hline(aes(yintercept=1e-13), linetype="dashed", color="#332288") +
    geom_hline(aes(yintercept=1e-5), linetype="dotdash", color="gray40") +
    geom_vline(aes(xintercept=(tstarTimepoint-29)), linetype="dotted") +
    annotate("segment", x=inferred_tau-29, xend=inferred_tau-29,
             y=10^-11, yend=10^-13, arrow = arrow(length = unit(.15,"cm"))) +
    annotate("segment", x=tstarTimepoint+50-29, xend=tstarTimepoint-29,
             y=10^-11, yend=10^-13, arrow=arrow(length = unit(.15,"cm"))) +
    # Annotate the length of the colonization delay.
    # annotate("segment", x=34-29, xend=tstarTimepoint-29,
    #          y=K*5, yend=K*5, arrow=arrow(ends = "both", angle = 90, length = unit(.1,"cm"))) +
    annotate("segment", x=lastUndetectableTimepoint-29, xend=tstarTimepoint-29,
             y=K*50, yend=K*50, arrow=arrow(ends = "both", angle = 90, length = unit(.1,"cm"))) +
    # Label the trajectory extrapolation.
    annotate("segment", x=lastUndetectableTimepoint-29, xend=inferred_tau-29,
             y=1e-5, yend=1e-13, linetype="dotted", 
             color=PALETTESTRAINS[unique(dataSpecies$typeOfRecoveryColonizationTurnover)]) +
    annotate(geom="rect", xmin=0, xmax=5, ymin=1e-13, ymax=1, fill=abxColor) +
    geom_line(aes(x=(timepoint), y=relative_abundance, group=(timePeriod),
                  color=factor(typeOfRecoveryColonizationTurnover))) +
    geom_point(aes(x=(timepoint), y=relative_abundance,
                   color=factor(typeOfRecoveryColonizationTurnover),
                   shape=factor(relative_abundance==1e-5)),
               alpha=0.8, size=1) +
    scale_color_manual(values=PALETTESTRAINS, name="Subject") +
    scale_shape_manual(values=c(16, 1)) +
    scale_y_continuous(trans='log10', limits=c(1e-13,5),
                       breaks=trans_breaks('log10', function(x) 10^x),
                       labels=trans_format('log10', math_format(10^.x))) +
    xlab("Study day") + ylab("Relative\nabundance") +
    DEFAULTS.THEME_PRINT +
    theme(plot.title=element_text(size=5, margin=margin(0,0,1,0)),
          strip.text.x=element_blank()) +
    guides(color = guide_legend(override.aes = list(size = 0.5), order=1),
           shape="none") +
    theme(legend.key.size = unit(0.5, "lines"),
          axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=5),
          plot.title=element_text(size=6, margin=margin(0,0,1,0)),) +
    ggtitle(paste0("Subject ", isubject, "\n", shortenSpeciesName(ispecies_id)))
  p4_species
  savePNGPDF(paste0(OUTDIR, "4-", isubject, "-", ispecies_id, "-main"),
             p4_species + guides(color="none"), height, width)
}

# Plot example colonization trajectories.
plotSpeciesTrajectory("XAA", "Alistipes_shahii_62199", 1.25, 2.3)
plotSpeciesTrajectory("XAA", "Alistipes_onderdonkii_55464", 1.25, 2.3)
plotSpeciesTrajectory("XDA", "Parabacteroides_johnsonii_55217", 1.25, 2.3)
plotSpeciesTrajectory("XKA", "Burkholderiales_bacterium_56577", 1.25, 2.3)


# Identify example colonization events with long colonization delays and short transit times.
dataTrajectoryFits <- read.table("workflow/analysis/fitColonizationTrajectories/out/speciesTrajectoriesFit-summary.txt",
                                 header=TRUE, stringsAsFactors = FALSE, sep="\t")
dataTrajectoryFitsOrdered <- dataTrajectoryFits %>%
  filter(tstarTimepoint>75) %>%
  arrange(timeIntervalToReachK, desc(timeOfRecoveryColonizationTurnover)) %>%
  dplyr::select(subject, species_id, timeIntervalToReachK, timeOfRecoveryColonizationTurnover)

# Import the species trajectory fits.
dataTrajectoryFits <- read.table("workflow/analysis/fitColonizationTrajectories/out/speciesTrajectoriesFit-summary.txt",
                                 header=TRUE, stringsAsFactors = FALSE, sep="\t")
# Summarize the events that occur during the follow-up study in each group of antibiotic-taking subjects.
# Exclude events that involve unknown strains from this visualization.
# In XBA, Anaerostipes_hadrus_55206 has a transit time of 15 days and saturation time of 749 days.
dataTrajectoryFitsFollowup <- dataTrajectoryFits %>%
  mutate(subjectResponse=annotateSubjectResponse(subject)) %>%
  filter(tstarTimepoint>75, typeOfRecoveryColonizationTurnover!="unknown strain",
         subjectResponse=="lasting response") %>%
  mutate(speciesColonizedTransientlyAll=
           ifelse(is.na(speciesColonizedTransientlyAll), FALSE, speciesColonizedTransientlyAll)) %>%
  group_by(subject, tstarTimepoint) %>%
  arrange(subject, tstarTimepoint,
          desc(typeOfRecoveryColonizationTurnover), (speciesColonizedTransientlyAll)) %>%
  mutate(timeDetectedResponseIndex=row_number())  %>%
  mutate(typeOfRecoveryColonizationTurnover=
           ifelse(speciesColonizedTransientlyAll, "new species, transient",
                  ifelse(typeOfRecoveryColonizationTurnover=="new species",
                         "new species, persistent", typeOfRecoveryColonizationTurnover)))
dataTrajectoryFitsFollowupDisplay <- dataTrajectoryFitsFollowup %>%
  mutate(tstarTimepoint=tstarTimepoint-29) %>%
  mutate(subjectResponseDisplay=annotateSubjectResponseNumSubjects(subject))

dataTrajectoryFitsFollowupDisplayTimepoints <-
  samplesRaw %>% filter(annotateSubjectResponse(subject)=="lasting response",
                        timepoint>75, !(sample %in% samplesXBAextraFollowup)) %>%
  dplyr::select(subject, timepoint) %>% mutate(timepoint=timepoint-29) %>% arrange(subject, timepoint)
dataTrajectoryFitsFollowupTimepoints <-
  as.character(sort(unique(c(dataTrajectoryFitsFollowupDisplayTimepoints$timepoint))))

# Plot a barplot of species recovery and colonization events in each subject response class.
p4_timelineColonizationBar <- dataTrajectoryFitsFollowupDisplay %>%
  #filter(subject!="XBA") %>%
  ggplot() +
  geom_blank(data=dataTrajectoryFitsFollowupDisplayTimepoints,
             aes(x=formatC(timepoint, width=3, format="d", flag="0"), y=0)) +
  geom_bar(aes(x=formatC(tstarTimepoint, width=3, format="d", flag="0"),
               fill=fct_relevel(typeOfRecoveryColonizationTurnover, rev(names(PALETTESTRAINS))),
               color=fct_relevel(typeOfRecoveryColonizationTurnover, rev(names(PALETTESTRAINS))))) +
  facet_grid(.~fct_relevel(subject, c("XKA","XEA","XDA","XAA","XBA")), scales="free", space="free_x") +
  scale_fill_manual(values=PALETTESTRAINSFILL, name="") +
  scale_color_manual(values=PALETTESTRAINS, name="") +
  scale_alpha_identity() + scale_y_continuous(breaks=pretty_breaks()) +
  theme(legend.key.size = unit(0.5, "lines")) + ylim(0,6) +
  xlab("Colonization/recovery time (days post-antibiotics)") + ylab("Number of\ncolonizers") +
  guides(fill = guide_legend(override.aes = list(size = 0.5)), color="none") +
  DEFAULTS.THEME_PRINT + theme(legend.position="bottom") +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=4.5),
        strip.background=element_rect(fill=paste0(PALETTE.ABXRESPONSE["lasting response"], "CC"))) +
  theme(legend.box.margin=margin(0,0,0,-10))
p4_timelineColonizationBar
savePNGPDF(paste0(OUTDIR, "4-lastingResponse-colonizers-all"),
           p4_timelineColonizationBar + guides(fill="none", color="none"), 1.25, 6.9)
savePNGPDF(paste0(OUTDIR, "4-lastingResponse-colonizers-legend"),
           get_legend(p4_timelineColonizationBar), 0.25, 3)


# Import the trajectory fit summary data.
dataSpeciesTrajectoriesFitSummary <- 
  read.table("workflow/analysis/fitColonizationTrajectories/out/speciesTrajectoriesFit-summary.txt",
             header=TRUE, stringsAsFactors = FALSE, sep="\t")
# Plot the distribution of the number of time intervals needed to go 
# from undetectable relative abundances to carrying capacity.
p4_transitTime_histogram <- dataSpeciesTrajectoriesFitSummary %>%
  filter(tstarTimepoint>70, typeOfRecoveryColonizationTurnover!="same strain",
         typeOfRecoveryColonizationTurnover!="unknown strain") %>%
  mutate(numTimepointsToReachKcondensed=ifelse(numTimepointsToReachK>=3,">2",as.character(numTimepointsToReachK))) %>%
  mutate(typeOfRecoveryColonizationTurnover=
           ifelse(!is.na(speciesColonizedTransientlyAll) & speciesColonizedTransientlyAll, "new species, transient",
                  ifelse(typeOfRecoveryColonizationTurnover=="new species",
                         "new species, persistent", typeOfRecoveryColonizationTurnover))) %>%
  group_by(typeOfRecoveryColonizationTurnover) %>%
  dplyr::count(numTimepointsToReachKcondensed) %>%
  ggplot() +
  geom_bar(aes(x=fct_relevel(numTimepointsToReachKcondensed, c("1","2",">2")), y=n,
               fill=fct_relevel(typeOfRecoveryColonizationTurnover, rev(names(PALETTESTRAINS))),
               color=fct_relevel(typeOfRecoveryColonizationTurnover, rev(names(PALETTESTRAINS)))), 
           stat="identity") +
  scale_fill_manual(values=PALETTESTRAINSFILL, name="", 
                    labels=c("new species,\ntransient", "new species,\npersistent", "new strain")) +
  scale_color_manual(values=PALETTESTRAINS, name="") +
  xlab("Transit time\n(number of timepoints)") + ylab("Number of\ncolonizers") +
  theme(legend.key.size = unit(0.5, "lines")) +
  guides(fill = guide_legend(override.aes = list(size = 0.5)), color="none") +
  DEFAULTS.THEME_PRINT + theme(legend.box.margin=margin(0,0,0,-10))
p4_transitTime_histogram
savePNGPDF(paste0(OUTDIR, "4-transitTime-histogram"), p4_transitTime_histogram,
           fig4height, 1.8)
# Calculate the proportion of colonization events in follow-up sampling that
# go from undetectability to carrying capacity in a single timepoint.
# Exclude events that involve the same strain and that involve unknown strains.
# 1 timepoint: 64/89 (72%); 2 timepoints: 14/89, 16%
dataSpeciesTrajectoriesFitSummary %>%
  filter(tstarTimepoint>70, typeOfRecoveryColonizationTurnover!="same strain",
         typeOfRecoveryColonizationTurnover!="unknown strain") %>%
  dplyr::count(numTimepointsToReachK) %>% mutate(totalEvents=sum(n), pctEvents=n/totalEvents)
# For colonization events with a transit time of a single timepoint,
# calculate the duration of that timepoint.
# min: 6 days, median: 80 days, max: 278 days
dataSpeciesTrajectoriesFitSummary %>%
  filter(tstarTimepoint>70, typeOfRecoveryColonizationTurnover!="same strain",
         typeOfRecoveryColonizationTurnover!="unknown strain",
           numTimepointsToReachK==1) %>%
  summarize(minTransitTime=min(timeIntervalToReachK), medianTransitTime=median(timeIntervalToReachK),
            maxTransitTime=max(timeIntervalToReachK))
# Calculate the distribution of the ratio of transit time / colonization delay
# and calculate the proportion of colonization events for which transit time
# is <25% of the observed colonization delay.
# For 45% of colonization events, transit time <25% of colonization delay.
dataSpeciesTrajectoriesFitSummary %>%
  filter(tstarTimepoint>70, typeOfRecoveryColonizationTurnover!="same strain",
         typeOfRecoveryColonizationTurnover!="unknown strain",
         numTimepointsToReachK==1) %>%
  mutate(transitTimePctColonizationDelay=timeIntervalToReachK/(tstarTimepoint-34)) %>%
  summarize(numOfColonizationDelay25=sum(transitTimePctColonizationDelay<0.25),
            totalEvents=n(), pctColonizationDelay25=numOfColonizationDelay25/totalEvents)
# Calculate the median inferred colonization time of colonization events
# detected during follow-up sampling.
# number of events with inferred tau > 6 months after end of abx: 44/112, 39%
dataSpeciesTrajectoriesFitSummary %>%
  filter(tstarTimepoint>70, typeOfRecoveryColonizationTurnover!="same strain") %>%
  summarize(numTauGreater6Mos=sum(inferred_tau>64), numTotalEvents=n(),
            pct=numTauGreater6Mos/numTotalEvents)


# Export the list of simultaneous colonization events in XDA.
# day 706, 5 events, previous timepoint 589, all absent at previous timepoint
# species: "Alistipes_putredinis_61533"       "Bacteroides_caccae_53434"         "Bacteroides_clarus_62282"
# "Parabacteroides_distasonis_56985" "Ruminococcus_bicirculans_59300"
dataSpeciesTrajectoriesFitSummary %>%
  filter(subject=="XDA", tstarTimepoint>70, typeOfRecoveryColonizationTurnover %in% c("new strain", "new species")) %>%
  group_by(tstarTimepoint) %>%
  mutate(numSimultaneousSpecies=n()) %>% ungroup() %>% filter(numSimultaneousSpecies==max(numSimultaneousSpecies)) %>%
  dplyr::select(species_id, tstarTimepoint, lastTimepointBeforeTstar, timeIntervalToReachK) %>%
  pull(species_id)
# Export the list of simultaneous colonization events in XAA.
# day 623, 7 events, previous timepoint 590, 6 absent at previous timepoint
# species: "Alistipes_onderdonkii_55464"            "Alistipes_putredinis_61533"
# "Alistipes_shahii_62199"                 "Bilophila_wadsworthia_57364"
# "Clostridiales_bacterium_61057"          "Faecalibacterium_cf_62236"
# "Intestinimonas_butyriciproducens_60001"
dataSpeciesTrajectoriesFitSummary %>%
  filter(subject=="XAA", tstarTimepoint>70) %>% group_by(tstarTimepoint) %>%
  mutate(numSimultaneousSpecies=n(), typeOfRecoveryColonizationTurnover %in% c("new strain", "new species")) %>%
  ungroup() %>% filter(numSimultaneousSpecies==max(numSimultaneousSpecies)) %>%
  dplyr::select(species_id, tstarTimepoint, lastTimepointBeforeTstar, timeIntervalToReachK) %>%
  pull(species_id)
# Export the list of simultaneous colonization events in XBA.
# day 388, 3 events, previous timepoint 380; 2 absent at previous timepoint
# species: "Dorea_longicatena_61473" "Erysipelotrichaceae_bacterium_59516" "Subdoligranulum_sp_62068"
dataSpeciesTrajectoriesFitSummary %>%
  filter(subject=="XBA", tstarTimepoint>70) %>% group_by(tstarTimepoint) %>%
  mutate(numSimultaneousSpecies=n()) %>% ungroup() %>% filter(numSimultaneousSpecies==max(numSimultaneousSpecies)) %>%
  dplyr::select(species_id, tstarTimepoint, lastTimepointBeforeTstar, timeIntervalToReachK) %>%
  pull(species_id)


# Count the number of shared species gained and lost in antibiotic-taking subjects.
dataSharedSpeciesGainedLost <-
  read.table("workflow/analysis/fitColonizationTrajectories/out/sharedSpecies-gainedLost.txt",
             header=TRUE, stringsAsFactors = FALSE)
p4_repeatability <- dataSharedSpeciesGainedLost %>%
  mutate(subjectResponse=annotateSubjectResponse(subject),
         subjectDisplay=ifelse(subjectResponse=="lasting response", subject, "other\nsubjects")) %>%
  pivot_longer(cols=starts_with("sharedSpecies"), names_to="speciesType", values_to="numSpecies") %>%
  mutate(speciesTypeDisplay=ifelse(speciesType=="sharedSpeciesLost","lost", "gained")) %>%
  #filter(subjectDisplay!="XBA") %>%
  ggplot() +
  geom_bar(aes(x=fct_relevel(subjectDisplay, c("XAA", "XBA", "XDA", "XEA", "XKA", "other\nsubjects")),
               y=numSpecies, fill=fct_rev(speciesTypeDisplay)), stat="identity",
           position=position_dodge(), width=0.7) +
  scale_fill_manual(values=c("#DDCC77","#117733"),
                    name="", labels=c("strains shared\nbefore abx,\nlost after abx",
                                      "strains not shared\nbefore abx,\ntransmitted after abx")) +
  xlab("Subject") + ylab("# of strains") +
  guides(fill = guide_legend(override.aes = list(size = 1))) +
  theme(legend.key.size = unit(0.5, "lines")) +
  DEFAULTS.THEME_PRINT +
  theme(axis.text.x=element_text(size=4.5)) +
  theme(legend.box.margin=margin(0,0,0,-10))
savePNGPDF(paste0(OUTDIR, "4-sharedStrainsGainedLost"), p4_repeatability, fig4height, 2.45)

# Calculate the relative abundance of Eubacterium eligens in XDA and XDB.
# This shared strain is lost after antibiotics.
# XDA: median pre-abx abundance: 2%
# XDB: median pre-abx abundance: 0.8%, median abundance: 2%
dataSpeciesTrajectories %>%
  filter(species_id=="Eubacterium_eligens_61678", hh=="XD") %>%
  group_by(subject) %>%
  summarize(medianPreAbxAbundance=median(relative_abundance[timepoint<29]),
            medianAbundance=median(relative_abundance),
            medianPostAbxAbundance=median(relative_abundance[timepoint>34]),
            maxPostAbxAbundance=max(relative_abundance[timepoint>34]))


# Figure 5 - Reemergence of pre-antibiotic strains ------------------------

fig5height <- 1.25
fig5mainwidth <- 1.5
fig5followupwidth <- 2
fig5followupwidth_narrow <- 1


source("workflow/analysis/plotHelpers.R")
OUTDIR <- "workflow/analysis/figures-v6/out/"

# Plot a large, labeled version of the species dynamics in XBA.
# Plot the community composition of XBA and XBB.
source("workflow/analysis/generateSpeciesAbundances/loadSpeciesAbundancePlots.R")
p5_XBAcomposition <- plotSpeciesAbundanceDaysFromAbxStart(speciesAbundances %>% group_by(species_id) %>%
                                            filter(subject=="XBA", max(relative_abundance)>LIMITOFDETECTION,
                                                   !(sample %in% samplesXBAextraFollowup)) %>%
                                            mutate(subject=annotateSubjectResponse(subject))) +
  DEFAULTS.THEME_PRINT
savePNGPDF(paste0(OUTDIR, "5-XBA-main"),
           p5_XBAcomposition + xlim(-30,35) + ylab("Relative abundance"), 1.5, 2.5)
savePNGPDF(paste0(OUTDIR, "5-XBA-followup"),
           p5_XBAcomposition + xlim(35,NA) + DEFAULTS.THEME_NOYAXIS, 1.5, 2.8)


plotStrainDynamics <- function(iData, iSubject, iPalette, axisIncrement, axisLabelIncrement) {
  annotationColor <- ifelse(iSubject=="XBA", "black", NA)
  maxRelativeAbundance <- max((iData %>% filter(subject==iSubject))$relative_abundance)
  yAxisLim <- ceiling(maxRelativeAbundance*(1/axisIncrement))/(1/axisIncrement)
  # Set axis breaks at every specified increment.
  axisBreaks <- seq(0,yAxisLim,axisLabelIncrement)
  # If there are more than five axis breaks, then subsample them.
  if(length(axisBreaks)>5){
    axisBreaks <- axisBreaks[seq(1, length(axisBreaks), 2)]
  }
  iData %>%
    filter(subject==iSubject) %>%
    ggplot() +
    geom_rect(aes(xmin=0, xmax=5, ymin=0, ymax=1.1*yAxisLim),
              fill=ifelse(iSubject=="XBA",abxColor,NA), alpha=0.8) +
    annotate(geom="point", x=436-29, y=maxRelativeAbundance,
             shape=25, fill=annotationColor, color=annotationColor) +
    geom_area(aes(x=timepoint-29, y=strainRelativeAbundance, fill=factor(strain))) +
    facet_wrap(~subject, ncol=1, scales="free") +
    scale_fill_manual(values=iPalette, name="Strain") +
    scale_y_continuous(limits=c(0, 1.1*yAxisLim),
                       breaks=axisBreaks) +
    guides(fill="none") +
    DEFAULTS.THEME_PRINT + theme(strip.text.x=element_text(size=6)) +
    theme(axis.title=element_blank(), axis.text=element_text(size=4.5))
}

# Import the strain dynamics for B. uniformis.
dataXBBuniformis <-
  read.table("workflow/analysis/inferXBAstrainDynamicsPersistenceRecolonization/out/strainAbundances-Bacteroides_uniformis_57318.txt",
             header=TRUE, stringsAsFactors = FALSE)
# Generate separate plots for XBA and XBB so that the axes can be scaled independently.
BuniformisPalette <- c("#882255","#332288","gray80")
names(BuniformisPalette) <- c("strain1","strain2","strainUnknown")
savePNGPDF(paste0(OUTDIR, "5-strainDynamics-Buniformis-XBA"),
           plotStrainDynamics(dataXBBuniformis, "XBA", BuniformisPalette, 0.05, 0.1) +
             xlim(-30,35),
           0.6*fig5height, fig5height)
savePNGPDF(paste0(OUTDIR, "5-strainDynamics-Buniformis-XBA-followup"),
           plotStrainDynamics(dataXBBuniformis, "XBA", BuniformisPalette, 0.05, 0.1) +
             xlim(35,NA) + DEFAULTS.THEME_NOYAXIS,
           0.6*fig5height, fig5height)
savePNGPDF(paste0(OUTDIR, "5-strainDynamics-Buniformis-XBB"),
           plotStrainDynamics(dataXBBuniformis, "XBB", BuniformisPalette, 0.05, 0.1) +
             xlim(-30,35),
           0.6*fig5height, fig5height)
savePNGPDF(paste0(OUTDIR, "5-strainDynamics-Buniformis-XBB-followup"),
           plotStrainDynamics(dataXBBuniformis, "XBB", BuniformisPalette, 0.05, 0.1) +
             xlim(35,NA) + DEFAULTS.THEME_NOYAXIS,
           0.6*fig5height, fig5height)

# Import the strain dynamics for B. vulgatus
dataXBBvulgatus <-
  read.table("workflow/analysis/inferXBAstrainDynamicsPersistenceRecolonization/out/strainAbundances-Bacteroides_vulgatus_57955.txt",
             header=TRUE, stringsAsFactors = FALSE)
# Generate separate plots for XBA and XBB so that the axes can be scaled independently.
BvulgatusPalette <- c("#882255","#CC6677","#44AA99","#882255","gray80")
names(BvulgatusPalette) <- c("strain1","strain2","strain3","strain4","strainUnknown")
savePNGPDF(paste0(OUTDIR, "5-strainDynamics-Bvulgatus-XBA"),
           plotStrainDynamics(dataXBBvulgatus, "XBA", BvulgatusPalette, 0.05, 0.2) +
             xlim(-30,35),
           0.6*fig5height, fig5height)
savePNGPDF(paste0(OUTDIR, "5-strainDynamics-Bvulgatus-XBA-followup"),
           plotStrainDynamics(dataXBBvulgatus, "XBA", BvulgatusPalette, 0.05, 0.2) +
             xlim(35,NA) + DEFAULTS.THEME_NOYAXIS,
           0.6*fig5height, fig5height)
savePNGPDF(paste0(OUTDIR, "5-strainDynamics-Bvulgatus-XBB"),
           plotStrainDynamics(dataXBBvulgatus, "XBB", BvulgatusPalette, 0.05, 0.2) +
             xlim(-30,35),
           0.6*fig5height, fig5height)
savePNGPDF(paste0(OUTDIR, "5-strainDynamics-Bvulgatus-XBB-followup"),
           plotStrainDynamics(dataXBBvulgatus, "XBB", BvulgatusPalette, 0.05, 0.2) +
             xlim(35,NA) + DEFAULTS.THEME_NOYAXIS,
           0.6*fig5height, fig5height)

# Import the strain dynamics for B. longum.
dataXBBlongum <-
  read.table("workflow/analysis/inferXBAstrainDynamicsPersistenceRecolonization/out/strainAbundances-Bifidobacterium_longum_57796.txt",
             header=TRUE, stringsAsFactors = FALSE)
# Generate separate plots for XBA and XBB so that the axes can be scaled independently.
BlongumPalette <- c("#882255","#332288","gray80")
names(BlongumPalette) <- c("strain1","strain2","strainUnknown")
savePNGPDF(paste0(OUTDIR, "5-strainDynamics-Blongum-XBA"),
           plotStrainDynamics(dataXBBlongum, "XBA", BlongumPalette, 0.05, 0.1) +
             xlim(-30,35),
           0.6*fig5height, fig5height)
savePNGPDF(paste0(OUTDIR, "5-strainDynamics-Blongum-XBA-followup"),
           plotStrainDynamics(dataXBBlongum, "XBA", BlongumPalette, 0.05, 0.1) +
             xlim(35,NA) + DEFAULTS.THEME_NOYAXIS,
           0.6*fig5height, fig5height)
savePNGPDF(paste0(OUTDIR, "5-strainDynamics-Blongum-XBB"),
           plotStrainDynamics(dataXBBlongum, "XBB", BlongumPalette, 0.05, 0.1) +
             xlim(-30,35),
           0.6*fig5height, fig5height)
savePNGPDF(paste0(OUTDIR, "5-strainDynamics-Blongum-XBB-followup"),
           plotStrainDynamics(dataXBBlongum, "XBB", BlongumPalette, 0.05, 0.1) +
             xlim(35,NA) + DEFAULTS.THEME_NOYAXIS,
           0.6*fig5height, fig5height)

# Import data on T6SS dynamics.
dataT6SS <- read.table("workflow/analysis/scratch/230305-T6SS/out/T6SScoverage-cleaned.txt",
                       header=TRUE, stringsAsFactors = FALSE)
# Exclude blacklisted samples from further analysis.
dataT6SS <- dataT6SS %>%
  filter((sample %in% samplesRaw$sample), !(sample %in% sampleBlacklist))
# Calculate the number of reads mapping to each bp of the T6SS genes
# per billion reads sequenced.
dataT6SS <- dataT6SS %>%
  mutate(normalizedAbundance=numReads/totalReads/geneLength*1e9)
# Extract the data for subjects XBA and XBB.
dataT6SSXB <- dataT6SS %>%
  filter(grepl("XB", sample)) %>%
  mutate(subject=substr(sample,1,3), timepoint=as.integer(substr(sample,5,7)))
# Extract the data for E/I genes only to distinguish more easily between genotypes.
# Annotate E/I genes based on their GA numbers as well.
dataT6SSXBEI <- dataT6SSXB %>%
  filter(geneFunction %in% c("E","I")) %>%
  mutate(genotype=paste0(GAtype,"-",geneFunction,geneNumber),
         genotypeShort=paste0(GAtype,"-",geneNumber))
# Retain only genotypes that are detected at some timepoint in these subjects.
dataT6SSXBEI<- dataT6SSXBEI %>%
  group_by(genotype) %>%
  filter(max(normalizedAbundance)>1)
DEFAULTS.PALETTE.TOL.COLORBLINDSAFE <- c("#332288", "#44AA99", "#88CCEE", "#DDCC77",
                                         "#CC6677", "#AA4499", "#882255", "#117733")
T6SSpalette <- c("#332288","#88CCEE","#882255","#117733","#44AA99")
names(T6SSpalette) <- sort(unique((dataT6SSXBEI %>% group_by(genotype) %>%
                                     filter(max(normalizedAbundance)>50))$genotypeShort))
source("workflow/analysis/generateSpeciesAbundances/loadSpeciesAbundancesFiltered.R")
plotT6SS <- function(iSubject){
  dataT6SSXBEI %>%
    group_by(genotype) %>%
    filter(max(normalizedAbundance)>50, subject==iSubject, geneFunction=="E",
           !(sample %in% samplesXBAextraFollowup), totalReads>1000) %>%
    ggplot() +
    geom_rect(aes(xmin=0, xmax=5, ymin=0, ymax=max(dataT6SSXBEI$normalizedAbundance)),
              fill=ifelse(iSubject=="XBA",abxColor,NA), alpha=0.8) +
    geom_line(aes(x=timepoint-29, y=normalizedAbundance, group=genotype,
                  color=factor(genotypeShort))) +
    facet_wrap(~subject, scales="free") +
    scale_color_manual(values=T6SSpalette, name="Effector") +
    scale_linetype_manual(values=c("solid","dashed"), name="Gene\ntype") +
    xlab("") + ylab("") + ylim(0,max(dataT6SSXBEI$normalizedAbundance)) +
    guides(color = guide_legend(override.aes = list(size = 0.5)),
           linetype = guide_legend(override.aes = list(size = 0.5))) +
    theme(legend.key.size = unit(0.5, "lines")) +
    DEFAULTS.THEME_PRINT
}

savePNGPDF(paste0(OUTDIR, "5-T6SS-XBA"),
           plotT6SS("XBA") + guides(linetype="none", color="none") + xlim(-30,35),
           fig5height, fig5mainwidth)
savePNGPDF(paste0(OUTDIR, "5-T6SS-XBA-followup"),
           plotT6SS("XBA") + guides(linetype="none", color="none") + xlim(35,NA) +
             DEFAULTS.THEME_NOYAXIS,
           fig5height, fig5followupwidth-0.25)
savePNGPDF(paste0(OUTDIR, "5-T6SS-XBB"),
           plotT6SS("XBB") + guides(linetype="none", color="none") + xlim(-30,35),
           fig5height, fig5mainwidth)
savePNGPDF(paste0(OUTDIR, "5-T6SS-XBB-followup"),
           plotT6SS("XBB") + guides(linetype="none", color="none") + xlim(35,NA) +
             DEFAULTS.THEME_NOYAXIS,
           fig5height, fig5followupwidth-0.25)
savePNGPDF(paste0(OUTDIR, "5-T6SS-legend"),
           get_legend(plotT6SS("XBB")), 0.75, 0.5)
