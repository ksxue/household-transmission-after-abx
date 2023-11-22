# Plot supplemental figures.

# Supp Fig 1 - Household demographics -------------------------------------

source("workflow/analysis/plotHelpers.R")
OUTDIR <- "workflow/analysis/supp-figures-v2/out/"

suppfig1height <- 1.25

# Import participant metadata.
dataParticipantMetadata <- 
  fread("workflow/analysis/parseParticipantMetadata/out/participantQuestionnairesCleaned.txt",
             stringsAsFactors = FALSE, header=TRUE, data.table=FALSE) 
# Rename columns to be more easily parsed by R.
colnames(dataParticipantMetadata) <- c("Subject", "Age", "Gender" ,"Hispanic",
                    "Race1", "Race2", "Race3", "Vegetarian", "Vegan", "GlutenFree", #10
                    "OtherDietaryRestrictions", "DairyIntolerance", "WheatIntolerance", #13
                    "EggIntolerance", "SoybeanIntolerance", "ShellfishIntolerance", "OtherIntolerance",
                    "LastAntibioticYrs", "LastTravelYrs", "LastInternationalTravelYrs", #20 
                    "TimeCurrentAddressMos", "NumAdults", "NumChildren", "NumDogs", "NumCats", "NumOtherPets" , #26
                    "Housemate1", "Housemate1LengthCohabitationMos", "Housemate1LengthCohabitationMosPairAvg", #29
                    "Housemate1ShareBed", "Housemate1ShareBedroom", "Housemate1ShareBathroom", #32
                    "Housemate1ShareKitchen", "Housemate1ShareDiningRoom", "Housemate1ShareCommonSpace", #35
                    "Housemate1Romantic", "Housemate1TimeSpentPct", "Housemate1TimeSpentPctPairAvg", #38
                    "Housemate1MealsPerWeek", "Housemate1MealsPerWeekPairAvg", #40
                    "Housemate2", "Housemate2LengthCohabitationMos", "Housemate2LengthCohabitationMosPairAvg", #43
                    "Housemate2ShareBed", "Housemate2ShareBedroom", "Housemate2ShareBathroom", #46
                    "Housemate2ShareKitchen", "Housemate2ShareDiningRoom", "Housemate2ShareCommonSpace", #49
                    "Housemate2Romantic", "Housemate2TimeSpentPct", "Housemate2TimeSpentPctPairAvg", #52
                    "Housemate2MealsPerWeek", "Housemate2MealsPerWeekPairAvg", #54
                    "Housemate3", "Housemate3LengthCohabitationMos", "Housemate3LengthCohabitationMosPairAvg", #57
                    "Housemate3ShareBed", "Housemate3ShareBedroom", "Housemate3ShareBathroom", #60
                    "Housemate3ShareKitchen", "Housemate3ShareDiningRoom", "Housemate3ShareCommonSpace", #63
                    "Housemate3Romantic", "Housemate3TimeSpentPct", "Housemate3TimeSpentPctPairAvg", #66
                    "Housemate3MealsPerWeek", "Housemate3MealsPerWeekPairAvg") #68

# Retain only subjects in study arm X.
dataParticipantMetadata <- dataParticipantMetadata %>%
  filter(Subject %in% subjectsX)

# Retain only demographic categories.
dataDemographics <- dataParticipantMetadata %>%
  dplyr::select(Subject, Gender, Hispanic, Race1, Race2, Race3)

# To match NIH reporting categories, create a multiracial category.
dataDemographics <- dataDemographics %>%
  mutate(Multiracial=ifelse(is.na(Race2),FALSE,TRUE)) %>%
  mutate(NIHRace=ifelse(Multiracial,"MoreThanOneRace",
                        ifelse(is.na(Race1),"NotReported",Race1))) %>%
  dplyr::select(Subject, Gender, Hispanic, NIHRace)

# To match NIH reporting categories, parse ethnicity.
dataDemographics <- dataDemographics %>%
  mutate(NIHEthnicity=ifelse(is.na(Hispanic),"NotReported",
                             ifelse(Hispanic=="None","NotHispanicLatino","HispanicLatino")))

# To match NIH reporting categories, parse self-reported gender.
dataDemographics <- dataDemographics %>%
  mutate(NIHGender=ifelse(!(Gender %in% c("Male","Female")),"NotReported", Gender))

# Plot participant genders.
# Group non-binary and not reported.
p1_gender <- dataDemographics %>%
  mutate(Gender=ifelse(is.na(Gender) | Gender=="Nonbinary","Not\nreported",Gender)) %>%
  ggplot() +
  geom_bar(aes(x=Gender)) +
  xlab("Gender") + ylab("Number of\nparticipants") +
  DEFAULTS.THEME_PRINT
savePNGPDF(paste0(OUTDIR, "1-gender"), p1_gender, suppfig1height, 1.5)

# Plot participant race and antibiotic responses.
p1_race <- dataDemographics %>%
  mutate(NIHRace=ifelse(NIHRace=="MoreThanOneRace","Multiracial",NIHRace)) %>%
  ggplot() +
  geom_bar(aes(x=fct_infreq(NIHRace))) +
  xlab("Race") + ylab("Number of\nparticipants") +
  DEFAULTS.THEME_PRINT
savePNGPDF(paste0(OUTDIR, "1-race"), p1_race, suppfig1height, 1.5)

# Plot participant ethnicity and antibiotic responses.
p1_ethnicity <- dataDemographics %>%
  mutate(NIHEthnicity= case_when(
    NIHEthnicity=="NotHispanicLatino" ~ "Not Hispanic\nor Latino",
    NIHEthnicity=="HispanicLatino" ~ "Hispanic\nor Latino",
    NIHEthnicity=="NotReported" ~ "Not\nreported"
  )) %>%
  ggplot() +
  geom_bar(aes(x=fct_infreq(NIHEthnicity))) +
  xlab("Ethnicity") + ylab("Number of\nparticipants") +
  DEFAULTS.THEME_PRINT +
  theme(axis.text.x=element_text(size=4))
savePNGPDF(paste0(OUTDIR, "1-ethnicity"), p1_ethnicity, suppfig1height, 1.5)

# Plot participant ages.
p1_age <- dataParticipantMetadata %>%
  ggplot() +
  geom_histogram(aes(x=Age), binwidth=5) +
  xlab("Age") + ylab("Number of\nparticipants") +
  DEFAULTS.THEME_PRINT
savePNGPDF(paste0(OUTDIR, "1-age"), p1_age, suppfig1height, 2)

# Import sample inventory.
# This lists all collected stool samples, regardless of whether or not they were sequenced.
sampleInventory <- read.table("data/householdCohort/sampleInventory-stool.txt",
                              header=TRUE, stringsAsFactors = FALSE) %>%
  filter(subject %in% subjectsX)
# Annotate samples that were sequenced.
sampleInventory <- sampleInventory %>%
  mutate(sequenced=sample %in% samplesRaw$sample,
         keyTimepoint=sample %in% samplesKeyTimepoints,
         sequencingDepth=ifelse(!sequenced, "0 Gbp",
                         ifelse(keyTimepoint, "10 Gbp", "2 Gbp")))
# Remove blacklisted samples.
sampleInventory <- sampleInventory %>%
  filter(!(sample %in% sampleBlacklist))
# Remove all follow-up samples that were not sequenced
# unless they were from subjects XBA and XBB,
# since the other follow-up samples had insufficient sequencing depth.
sampleInventory <- sampleInventory %>%
  filter(timepoint<75 | sample %in% samplesRaw$sample |
           substr(subject,1,2)=="XB")

# Generate a plot that shows the samples collected and sequenced,
# as well as their target sequencing depth.
PALETTE.SEQDEPTH <- c("gray30","#CC6677","#882255")
names(PALETTE.SEQDEPTH) <- c("0 Gbp", "2 Gbp", "10 Gbp")
p1_inventoryMain <- sampleInventory %>% 
  group_by(subject) %>%
  mutate(timepoint=timepoint-29) %>%
  mutate(minTimepointMain=min(timepoint),
         maxTimepointMain=max(timepoint[timepoint<45]),
         timePeriod=ifelse(timepoint<45, "main study", "follow-up sampling")) %>% 
  filter(timePeriod=="main study") %>%
  arrange(fct_relevel(sequencingDepth, c("0 Gbp","2 Gbp","10 Gbp")), subject, timepoint) %>%
  ggplot() +
  geom_segment(aes(y=fct_rev(subject), yend=fct_rev(subject), 
                   x=minTimepointMain, xend=maxTimepointMain),
               linewidth=0.4) +
  geom_point(aes(x=timepoint, y=fct_rev(subject), 
                 color=fct_relevel(sequencingDepth, c("0 Gbp","2 Gbp","10 Gbp")),
                 size=fct_relevel(sequencingDepth, c("0 Gbp","2 Gbp","10 Gbp")))) +
  facet_wrap(~timePeriod) +
  scale_color_manual(values=PALETTE.SEQDEPTH, name="Sequencing\ndepth") +
  scale_size_manual(values=c(0.5,1,1.5), name="Sequencing\ndepth") +
  xlab("Study day") + ylab("Subject") +
  DEFAULTS.THEME_PRINT +
  guides(color = guide_legend(override.aes = list(size = 1)), size="none") +
  theme(legend.key.size = unit(0.5, "lines"),
        axis.text.y=element_text(size=5))
savePNGPDF(paste0(OUTDIR, "1-inventory-main"), 
           p1_inventoryMain + guides(color="none"), 4.5, 3)

p1_inventoryFollowup <- sampleInventory %>% 
  group_by(subject) %>%
  mutate(timepoint=timepoint-29) %>%
  mutate(minTimepointFollowup=min(timepoint[timepoint>45]),
         maxTimepointFollowup=max(timepoint),
         timePeriod=ifelse(timepoint<45, "main study", "follow-up sampling")) %>%
  filter(timePeriod=="follow-up sampling") %>%
  arrange(fct_relevel(sequencingDepth, c("0 Gbp","2 Gbp","10 Gbp")), subject, timepoint) %>%
  ggplot() +
  geom_segment(aes(y=fct_rev(subject), yend=fct_rev(subject), 
                   x=minTimepointFollowup, xend=maxTimepointFollowup),
               linewidth=0.4) +
  geom_point(aes(x=timepoint, y=fct_rev(subject), 
                 color=fct_relevel(sequencingDepth, c("0 Gbp","2 Gbp","10 Gbp")),
                 size=fct_relevel(sequencingDepth, c("0 Gbp","2 Gbp","10 Gbp")))) +
  facet_wrap(~timePeriod) +
  scale_color_manual(values=PALETTE.SEQDEPTH, name="Sequencing\ndepth") +
  scale_size_manual(values=c(0.5,1,1.5), name="Sequencing\ndepth") +
  xlim(45,NA) +
  xlab("Study day") + ylab("Subject") +
  DEFAULTS.THEME_PRINT +
  guides(color = guide_legend(override.aes = list(size = 1)), size="none") +
  theme(legend.key.size = unit(0.5, "lines"),
        axis.text.y=element_text(size=5),
        axis.title.y=element_blank())
savePNGPDF(paste0(OUTDIR, "1-inventory-followup"), 
           p1_inventoryFollowup + guides(color="none"), 4.5, 3)
savePNGPDF(paste0(OUTDIR, "1-inventory-legend"), 
           get_legend(p1_inventoryMain), 0.75, 0.75)


# Supp figure 2 - identifying shared strains ------------------------------

source("workflow/analysis/plotHelpers.R")
OUTDIR <- "workflow/analysis/supp-figures-v2/out/"

suppfig2height <- 1.25

# Import fixed differences data on the similarity thresholds for each species.
dataFixedDiffsThresholds <- read.table("workflow/analysis/identifySharedStrains-fixedDiffs/out/fixedDiffs-btwnHouseholdThresholds.txt",
                                       header=TRUE, stringsAsFactors = FALSE)
p2_fixedDiffsThresholds <- dataFixedDiffsThresholds %>%
  filter(!is.na(pctGenomeLength)) %>%
  ggplot() +
  #geom_vline(aes(xintercept=median(pctGenomeLength, na.rm=TRUE)), linetype="dashed") +
  geom_histogram(aes(x=pctGenomeLength), bins=15) +
  scale_x_log10(limits=c(1e-6,1)) +
  xlab("% genome length") + ylab("Number of species") +
  DEFAULTS.THEME_PRINT
savePNGPDF(paste0(OUTDIR, "2-fixedDiffsThresholds"), p2_fixedDiffsThresholds,
           suppfig2height, 1.75)

# Import strain fishing data for sample species.
dataStrainFishingSampleSpecies <- read.table("workflow/analysis/identifySharedStrains-strainFishing/out/strainFishing-sampleSpecies.txt.gz",
                                             header=TRUE, stringsAsFactors = FALSE)
# Plot the distribution of % private SNPs detected and the 75% threshold.
p2_strainFishingSampleSpecies <- dataStrainFishingSampleSpecies %>%
  mutate(speciesShort=shortenSpeciesName(species),
         typeDisplay=case_when(
           type=="sameSubject" ~ "same\nsubject",
           type=="sameHousehold" ~ "same\nhousehold",
           type=="diffHousehold" ~ "different\nhousehold"
         )) %>%
  ggplot() +
  geom_violin(aes(x=fct_rev(typeDisplay), y=pctSitesDetected, 
                  fill=fct_rev(typeDisplay), color=fct_rev(typeDisplay)), 
              scale="width") +
  geom_hline(aes(yintercept=0.75), linetype="dashed") +
  scale_fill_manual(name="Population\npair", values=c("#AA4499", "#44AA99", "#332288")) +
  scale_color_manual(name="Population\npair", values=c("#AA4499", "#44AA99", "#332288")) +
  facet_wrap(~speciesShort) +
  xlab("Population pair") + ylab("% of private SNPs\ndetected") +
  DEFAULTS.THEME_PRINT +
  guides(color = guide_legend(override.aes = list(size = 1))) +
  theme(legend.key.size = unit(0.5, "lines"), axis.text.x=element_text(size=4.5))
savePNGPDF(paste0(OUTDIR, "2-strainFishing-sampleSpecies-violin"), p2_strainFishingSampleSpecies,
           suppfig2height, 5)

# Import fixed differences data for sample species.
dataFixedDiffs <- read.table("workflow/analysis/identifySharedStrains-fixedDiffs/out/fixedDiffs-allSamples.txt.gz",
                             header=TRUE, stringsAsFactors = FALSE)
# Retain only data for sample species of interest.
dataFixedDiffs <- dataFixedDiffs %>%
  filter(species %in% c("Bacteroides_fragilis_54507", "Bifidobacterium_longum_57796",
                        "Eubacterium_eligens_61678"))
# Plot the distribution of fixed differences and the 1% FDR threshold for each species.
p2_fixedDiffsSampleSpecies <- dataFixedDiffs %>%
  mutate(speciesShort=shortenSpeciesName(species),
         typeDisplay=case_when(
           type=="sameSubject" ~ "same\nsubject",
           type=="sameHousehold" ~ "same\nhousehold",
           type=="diffHousehold" ~ "different\nhousehold"
         )) %>%
  ggplot() +
  geom_violin(aes(x=fct_rev(typeDisplay), y=fixedDiffs+1, 
                  fill=fct_rev(typeDisplay), color=fct_rev(typeDisplay)), 
              scale="width") +
  geom_hline(data=dataFixedDiffs %>% dplyr::select(species, btwnThreshold) %>% unique() %>%
               mutate(speciesShort=shortenSpeciesName(species)),
             aes(yintercept=btwnThreshold), linetype="dashed") +
  scale_fill_manual(name="Population\npair", values=c("#AA4499", "#44AA99", "#332288")) +
  scale_color_manual(name="Population\npair", values=c("#AA4499", "#44AA99", "#332288")) +
  facet_wrap(~speciesShort) +
  scale_y_log10() +
  xlab("Population pair") + ylab("Number of\nfixed differences") +
  DEFAULTS.THEME_PRINT +
  guides(color = guide_legend(override.aes = list(size = 1))) +
  theme(legend.key.size = unit(0.5, "lines"), axis.text.x=element_text(size=4.5))
savePNGPDF(paste0(OUTDIR, "2-fixedDiffs-sampleSpecies-violin"), p2_fixedDiffsSampleSpecies,
           suppfig2height, 5)

# Import strain fishing data on the proportion of comparisons that involve cousins.
dataStrainFishingPctCousins <- read.table("workflow/analysis/identifySharedStrains-strainFishing/out/strainFishing-pctCousins.txt",
                                          header=TRUE, stringsAsFactors = FALSE)
p2_strainFishingPctCousins <- dataStrainFishingPctCousins %>%
  ggplot() +
  geom_histogram(aes(x=pctCousins), bins=30) +
  xlim(NA, 0.3) +
  xlab("% of between-household pairs,\n>75% private SNPs detected") + 
  ylab("Number of species") +
  DEFAULTS.THEME_PRINT
savePNGPDF(paste0(OUTDIR, "2-strainFishingPctCousins"), p2_strainFishingPctCousins,
           suppfig2height, 1.75)
# 19% of species (16 of 84) have 0 detected cousins.
nrow(dataStrainFishingPctCousins %>% filter(pctCousins==0))/nrow(dataStrainFishingPctCousins)


# Import a comparison of the shared strain calls from fixed differences and strain sharing.
dataStrainSharingComparison <- 
  read.table("workflow/analysis/compareSharedStrainCalls/out/strainSharingCalls-consolidated.txt",
             header=TRUE, stringsAsFactors = FALSE)
# Filter to include only comparisons of population pairs from the initial sampling timepoints
# in subjects in study arm X.
dataStrainSharingComparisonInitial <- dataStrainSharingComparison %>%
  filter(sample1 %in% samplesInitial, sample2 %in% samplesInitial,
         subject1 %in% subjectsX, subject2 %in% subjectsX)
dataStrainSharingComparisonInitial <- dataStrainSharingComparisonInitial %>%
  mutate(typeConcordance=case_when(
    is.na(sharedFixedDiffs) ~ "strain fishing only",
    is.na(sharedStrainFishing) ~ "fixed diffs only",
    sharedFixedDiffs==sharedStrainFishing ~ "concordant",
    sharedFixedDiffs ~ "discordant,\nfixed diffs TRUE",
    sharedStrainFishing ~ "discordant,\nstrain fishing TRUE"
  ),
  type=case_when(
    is.na(sharedFixedDiffs) ~ "strain fishing only",
    is.na(sharedStrainFishing) ~ "fixed diffs only",
    !is.na(sharedFixedDiffs) & !is.na(sharedStrainFishing) ~ "both methods"
  ))
dataStrainSharingComparisonInitialSummary <- dataStrainSharingComparisonInitial %>%
  arrange(fct_relevel(type, c("fixed diffs only", "both methods", "strain fishing only"))) %>%
  group_by(type) %>% mutate(numPairs=n()) %>%
  dplyr::select(type, numPairs) %>% unique() %>% ungroup() %>%
  mutate(labelPosition=0.5*numPairs+cumsum(lag(numPairs, default=0)))
dataStrainSharingComparisonInitialConcordanceSummary <- dataStrainSharingComparisonInitial %>%
  arrange(fct_relevel(typeConcordance, 
                      c("fixed diffs only", "discordant,\nfixed diffs TRUE", "concordant", 
                        "discordant,\nstrain fishing TRUE", "strain fishing only"))) %>%
  group_by(typeConcordance) %>% mutate(numPairs=n()) %>%
  dplyr::select(typeConcordance, numPairs) %>% unique() %>% ungroup() %>%
  mutate(labelPosition=0.5*numPairs+cumsum(lag(numPairs, default=0)))
p2_callSets <- dataStrainSharingComparisonInitialSummary %>%
  ggplot() +
  geom_bar(aes(x="type", y=numPairs,
               fill=fct_relevel(type, rev(c("fixed diffs only", "both methods", "strain fishing only")))),
           color=NA, stat="identity") +
  geom_text(aes(x="type", y=labelPosition, 
                label=paste0(type,"\n",numPairs," pairs")), size=2, color="white") +
  scale_fill_manual(values=c("gray30","gray","gray30"), name="") +
  coord_flip() +
  DEFAULTS.THEME_PRINT +
  theme_nothing()
p2_callConcordance <- dataStrainSharingComparisonInitialConcordanceSummary %>%
  ggplot() +
  scale_x_discrete() +
  geom_bar(aes(x="type", y=numPairs,
               fill=fct_relevel(typeConcordance, 
                                rev(c("fixed diffs only", "discordant,\nfixed diffs TRUE", "concordant", 
                                  "discordant,\nstrain fishing TRUE", "strain fishing only")))),
           color=NA, stat="identity") +
  geom_text(data = dataStrainSharingComparisonInitialConcordanceSummary %>%
              filter(!grepl("only",typeConcordance)),
            aes(x=ifelse(grepl("discordant",typeConcordance), 1, 1), y=labelPosition, 
                label=paste0(typeConcordance,"\n",numPairs," pairs")), size=2, color="white") +
  scale_fill_manual(values=c("gray30","#882255","#44AA99","#882255","gray30"), name="") +
  coord_flip() +
  DEFAULTS.THEME_PRINT +
  theme_nothing()
p2_compareCallMethods_blank <- 
  (dataStrainSharingComparisonInitialConcordanceSummary %>%
     ggplot() +
     scale_x_discrete() +
     geom_bar(aes(x="type", y=numPairs,
                  fill=fct_relevel(typeConcordance, 
                                   rev(c("fixed diffs only", "discordant,\nfixed diffs TRUE", "concordant", 
                                         "discordant,\nstrain fishing TRUE", "strain fishing only")))),
              color=NA, stat="identity") +
     scale_fill_manual(values=c("gray30","#882255","#44AA99","#882255","gray30"), name="") +
     coord_flip() +
     DEFAULTS.THEME_PRINT +
     theme_nothing())
savePNGPDF(paste0(OUTDIR, "2-compareCallMethods"), p2_callSets / p2_callConcordance,
           suppfig2height, 7)
savePNGPDF(paste0(OUTDIR,"2-compareCallMethods-blank"),p2_compareCallMethods_blank,
           0.5, 3)

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
# Order subject pairs based on the relative abundance of shared strains.
samplePairOrder <- 
  dataRelAbundanceStrainSharingSummaryRaw %>%
  filter(sample %in% samplesInitial, annotationSample %in% samplesInitial,
         sample!=annotationSample, sample %in% samplesX, annotationSample %in% samplesX,
         substr(sample,1,3)!=substr(annotationSample,1,3)) %>%
  mutate(samplePair=paste0(sample,"_",annotationSample)) %>%
  complete(samplePair, annotation, fill=list(numSpecies=0, totRelAbundance=0)) %>%
  mutate(sample=substr(samplePair,1,7), annotationSample=substr(samplePair,9,16)) %>%
  filter(annotation=="strainsShared") %>%
  mutate(hh=substr(sample,1,2)) %>% group_by(hh) %>% 
  mutate(meanTotRelAbundance=mean(totRelAbundance)) %>%
  mutate(subjectPair=paste0(substr(sample,1,3),"-",substr(annotationSample,1,3))) %>%
  arrange(meanTotRelAbundance, totRelAbundance, samplePair)
# Plot the relative abundance of shared species and strains.
p2_stackedBarPlots <- dataRelAbundanceStrainSharing %>%
  ggplot() +
  geom_bar(aes(x=fct_relevel(samplePair, samplePairOrder$samplePair), y=relative_abundance, 
               fill=fct_relevel(annotation, c("speciesNotShared","strainsUnknown","strainsNotShared","strainsShared"))), 
           stat="identity", color="black", linewidth=0.1) +
  scale_fill_manual(values=sharingAnnotationPalette, name="Sharing status",
                    labels=c("species not shared", "strains unknown", "strains not shared", "strains shared")) +
  xlab("Subject pair") + ylab("Relative abundance") +
  DEFAULTS.THEME_PRINT +
  guides(fill = guide_legend(override.aes = list(size = 0.5))) +
  scale_x_discrete(labels=samplePairOrder$subjectPair) +
  theme(legend.key.size = unit(0.5, "lines")) +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=4)) +
  theme(axis.ticks.x=element_blank()) +
  theme(axis.line=element_blank()) +
  theme(axis.text.y=element_text(size=6))
p2_stackedBarPlots
savePNGPDF(paste0(OUTDIR, "2-sharingAnnotations-relAbundance-stackedBars"),
          p2_stackedBarPlots, 1.75, 6.8)

# Plot the proportion of strains shared in non-cohabiting subjects.
# Import strain sharing annotation summaries from the initial timepoint.
dataStrainSharingInitialSummaryRaw <- 
  read.table("workflow/analysis/integrateSpeciesAbundancesSharedStrains/out/annotatedSpeciesAbundances-initialTimepoints-summary.txt",
             header=TRUE, stringsAsFactors = FALSE)
# Use complete to infer when subject pairs have no strains in common.
dataStrainSharingInitialSummary <- dataStrainSharingInitialSummaryRaw %>%
  mutate(samplePair=paste0(sample,"_",annotationSample)) %>%
  complete(samplePair, annotation, fill = list(totRelAbundance=0, numSpecies=0)) %>%
  mutate(sample=substr(samplePair,1,7), annotationSample=substr(samplePair,9,16)) %>%
  filter(sample!=annotationSample) %>%
  mutate(hh=ifelse(substr(sample,1,2)==substr(annotationSample,1,2), "sameHh","diffHh"))

# Compare the amount of strain sharing for cohabiting and non-cohabiting pairs.
p_hhcomparisonviolin <- dataStrainSharingInitialSummary %>%
  mutate(displayAnnotation=displayAnnotations[annotation],
         displayHh=ifelse(hh=="sameHh","same household", "different household")) %>%
  ggplot() +
  geom_violin(aes(x=displayAnnotation, y=totRelAbundance,
                  fill=factor(annotation)), scale="width") +
  geom_boxplot(aes(x=displayAnnotation, y=totRelAbundance, group=factor(annotation)), 
               width=0.1, outlier.size=0.1) +
  facet_wrap(~fct_rev(displayHh)) +
  scale_fill_manual(values=sharingAnnotationPalette) +
  guides(fill="none") +
  xlab("") + ylab("Total\nrelative abundance") +
  DEFAULTS.THEME_PRINT + theme(axis.text.x=element_text(size=5))
savePNGPDF(paste0(OUTDIR, "2-withinBtwnComparison-violin"), p_hhcomparisonviolin, 1.25, 3.5)


# Supp figure 3 - abx responses -------------------------------------------

source("workflow/analysis/plotHelpers.R")
OUTDIR <- "workflow/analysis/supp-figures-v2/out/"

suppfig3height <- 1.25

# Import FST values calculated across key timepoints to compare subjects.
# Import FSTruct data.
dataFstructKeyTimepoints <- 
  read.table("workflow/analysis/calculateDiversityStats/out/speciesCV-keyTimepoints.txt",
             header=TRUE, stringsAsFactors = FALSE) %>%
  mutate(subjectTreatment=annotateSubjectTreatment(subject),
         subjectResponse=annotateSubjectResponse(subject)) 
# Plot a comparison of antibiotic-taking and control subjects.
p3_FST_main <- dataFstructKeyTimepoints %>%
  arrange(desc(subjectResponse)) %>%
  ggplot() +
  geom_point(aes(x=fct_rev(subjectTreatment), y=ratio, color=factor(subjectResponse)),
             position=position_jitter(height=0, width=0.2, seed=0), alpha=0.8, size=1) +
  scale_color_manual(values=PALETTE.ABXRESPONSE) +
  xlab("Subject") + ylab("Compositional\nvariability\n(normalized FST)") +
  guides(color="none") +
  DEFAULTS.THEME_PRINT
savePNGPDF(paste0(OUTDIR, "3-FST-keyTimepoints"), p3_FST_main, suppfig3height, suppfig3height)
# Statistically compare the values of FAVA in control and antibiotic-taking subjects.
# Wilcoxon rank sum test p=3.312e-6
wilcox.test(dataFstructKeyTimepoints %>% filter(subjectTreatment=="control") %>% pull(ratio),
            dataFstructKeyTimepoints %>% filter(subjectTreatment=="abx") %>% pull(ratio))

# Plot alpha diversity (Shannon effective species) at different timepoints.
dataAlphaEffectiveSpecies <- read.table("workflow/analysis/calculateDiversityStats/out/speciesAlpha.txt",
                                        header=TRUE, stringsAsFactors = FALSE) %>%
  filter(alphaStat=="ShannonEffectiveSpecies", type=="filtSample",
         sample %in% samplesXmain) %>%
  mutate(subject=substr(sample,1,3), timepoint=as.integer(substr(sample,5,7)),
         subjectResponse=annotateSubjectResponse(subject))
dataAlphaEffectiveSpeciesSummary <- dataAlphaEffectiveSpecies %>%
  group_by(subject, subjectResponse) %>% 
  summarize(medianPreAbxDiversity=median(value[timepoint<29]),
         minPostAbxDiversity=min(value[timepoint>=29 & timepoint<75]),
         maxTimepoint=max(timepoint[timepoint<75]),
         finalDiversityMainStudy=value[timepoint==maxTimepoint]) %>%
  mutate(subjectTreatment=annotateSubjectTreatment(subject))

p3_medianMinimumDiversity <- dataAlphaEffectiveSpeciesSummary %>%
  arrange(desc(subjectTreatment), desc(subjectResponse)) %>%
  ggplot() +
  geom_abline(aes(intercept=0, slope=1), linetype="dashed") +
  geom_point(aes(x=medianPreAbxDiversity, y=minPostAbxDiversity,
                 color=fct_relevel(subjectResponse, ABXRESPONSE)), size=0.5, alpha=0.8) +
  scale_color_manual(values=PALETTE.ABXRESPONSE, name="Subject\nresponse") +
  xlim(0,NA) + ylim(0,NA) +
  guides(color = guide_legend(override.aes = list(size = 0.5))) +
  theme(legend.key.size = unit(0.5, "lines"),
        legend.box.margin=margin(0,0,0,-10)) +
  xlab("Median pre-abx\ndiversity") + ylab("Minimum post-abx\ndiversity") +
  DEFAULTS.THEME_PRINT
p3_medianFinalDiversity <- dataAlphaEffectiveSpeciesSummary %>%
  arrange(desc(subjectTreatment), desc(subjectResponse)) %>%
  ggplot() +
  geom_abline(aes(intercept=0, slope=1), linetype="dashed") +
  geom_point(aes(x=medianPreAbxDiversity, y=finalDiversityMainStudy,
                 color=fct_relevel(subjectResponse, ABXRESPONSE)), size=0.5, alpha=0.8) +
  scale_color_manual(values=PALETTE.ABXRESPONSE, name="Subject\nresponse") +
  xlim(0,NA) + ylim(0,NA) +
  guides(color = guide_legend(override.aes = list(size = 0.5))) +
  theme(legend.key.size = unit(0.5, "lines"),
        legend.box.margin=margin(0,0,0,-10)) +
  xlab("Median pre-abx\ndiversity") + ylab("Final diversity,\n1 mo post-abx") +
  DEFAULTS.THEME_PRINT
savePNGPDF(paste0(OUTDIR, "3-alphaDiversity-medianMinimum"), 
           p3_medianMinimumDiversity + guides(color="none"), suppfig3height, suppfig3height+0.1)
savePNGPDF(paste0(OUTDIR, "3-alphaDiversity-medianFinal"), 
           p3_medianFinalDiversity + guides(color="none"), suppfig3height, suppfig3height+0.1)


# Plot family-level JSD from the initial timepoint.
# Import the family-level JSD data.
dataFamilyJSD <- read.table("workflow/analysis/calculateDiversityStats/out/familiesBeta.txt.gz",
                            header=TRUE, stringsAsFactors = FALSE)
# Include only the distances calculated from the initial timepoint.
dataFamilyJSDFromInitial <- dataFamilyJSD %>%
  filter(substr(sample1,1,3)==substr(sample2,1,3), method=="jsd",
         sample1 %in% samplesInitial) %>%
  mutate(timepoint=as.integer(substr(sample2,5,7)),
         subject=substr(sample1,1,3)) %>%
  filter(subject %in% subjectsX)
# Calculate the distribution of family-level JSD between samples from
# unrelated individuals at the initial timepoint.
dataFamilyJSDNonCohabiting <- dataFamilyJSD %>%
  filter(substr(sample1,1,2)!=substr(sample2,1,2), method=="jsd",
         sample1 %in% samplesInitial, sample2 %in% samplesInitial) %>%
  mutate(subject1=substr(sample1,1,3), subject2=substr(sample2,1,3)) %>%
  filter(subject1 %in% subjectsX, subject2 %in% subjectsX)
familyJSDNonCohabiting <- quantile(dataFamilyJSDNonCohabiting$value, c(0.025,0.5,0.975))

# Plot the maximum and final JSD for each subject.
dataFamilyJSDFromInitialMaxFinal <- dataFamilyJSDFromInitial %>%
  group_by(subject) %>%
  filter(timepoint<75) %>%
  summarize(maxJSD=max(value), finalJSD=value[timepoint==max(timepoint)])
p3_JSDmaxfinal <- dataFamilyJSDFromInitialMaxFinal %>%
  mutate(subjectResponse=annotateSubjectResponse(subject),
         subjectTreatment=annotateSubjectTreatment(subject)) %>%
  arrange(desc(subjectTreatment), desc(subjectResponse)) %>%
  ggplot() +
  geom_abline(intercept=0, slope=1, linetype="dashed") +
  geom_point(aes(x=maxJSD, y=finalJSD, 
                 color=fct_relevel(subjectResponse, ABXRESPONSE)), alpha=0.8, size=1) + 
  scale_color_manual(values=PALETTE.ABXRESPONSE,
                     name="Subject\nresponse") +
  xlim(0,0.7) + ylim(0,0.7) +
  guides(color = guide_legend(override.aes = list(size = 0.75))) +
  theme(legend.key.size = unit(0.5, "lines")) +
  xlab("Maximum\nfamily-level JSD") + ylab("Final\nfamily-level JSD") +
  DEFAULTS.THEME_PRINT
savePNGPDF(paste0(OUTDIR, "3-familyJSD-maxfinal"), p3_JSDmaxfinal, 1.25, 2.5)


# Import participant metadata.
dataParticipantMetadata <- 
  fread("workflow/analysis/parseParticipantMetadata/out/participantQuestionnairesCleaned.txt",
        stringsAsFactors = FALSE, header=TRUE, data.table=FALSE) 
# Rename columns to be more easily parsed by R.
colnames(dataParticipantMetadata) <- c("Subject", "Age", "Gender" ,"Hispanic",
                                       "Race1", "Race2", "Race3", "Vegetarian", "Vegan", "GlutenFree", #10
                                       "OtherDietaryRestrictions", "DairyIntolerance", "WheatIntolerance", #13
                                       "EggIntolerance", "SoybeanIntolerance", "ShellfishIntolerance", "OtherIntolerance",
                                       "LastAntibioticYrs", "LastTravelYrs", "LastInternationalTravelYrs", #20 
                                       "TimeCurrentAddressMos", "NumAdults", "NumChildren", "NumDogs", "NumCats", "NumOtherPets" , #26
                                       "Housemate1", "Housemate1LengthCohabitationMos", "Housemate1LengthCohabitationMosPairAvg", #29
                                       "Housemate1ShareBed", "Housemate1ShareBedroom", "Housemate1ShareBathroom", #32
                                       "Housemate1ShareKitchen", "Housemate1ShareDiningRoom", "Housemate1ShareCommonSpace", #35
                                       "Housemate1Romantic", "Housemate1TimeSpentPct", "Housemate1TimeSpentPctPairAvg", #38
                                       "Housemate1MealsPerWeek", "Housemate1MealsPerWeekPairAvg", #40
                                       "Housemate2", "Housemate2LengthCohabitationMos", "Housemate2LengthCohabitationMosPairAvg", #43
                                       "Housemate2ShareBed", "Housemate2ShareBedroom", "Housemate2ShareBathroom", #46
                                       "Housemate2ShareKitchen", "Housemate2ShareDiningRoom", "Housemate2ShareCommonSpace", #49
                                       "Housemate2Romantic", "Housemate2TimeSpentPct", "Housemate2TimeSpentPctPairAvg", #52
                                       "Housemate2MealsPerWeek", "Housemate2MealsPerWeekPairAvg", #54
                                       "Housemate3", "Housemate3LengthCohabitationMos", "Housemate3LengthCohabitationMosPairAvg", #57
                                       "Housemate3ShareBed", "Housemate3ShareBedroom", "Housemate3ShareBathroom", #60
                                       "Housemate3ShareKitchen", "Housemate3ShareDiningRoom", "Housemate3ShareCommonSpace", #63
                                       "Housemate3Romantic", "Housemate3TimeSpentPct", "Housemate3TimeSpentPctPairAvg", #66
                                       "Housemate3MealsPerWeek", "Housemate3MealsPerWeekPairAvg") #68

# Retain only subjects in study arm X.
dataParticipantMetadata <- dataParticipantMetadata %>%
  filter(Subject %in% subjectsX)

# Plot the distribution of prior antibiotic usage and antibiotic responses.
p3_recentAbx <- dataParticipantMetadata %>%
  mutate(LastAbxDisplay=case_when(
    is.na(LastAntibioticYrs) ~ "Unknown",
    LastAntibioticYrs==0.25 ~ "<6 months",
    LastAntibioticYrs==0.67 ~ "6-12 months",
    LastAntibioticYrs==1.50 ~ "1-2 years",
    LastAntibioticYrs==4.00 ~ "3-5 years",
    LastAntibioticYrs==5.00 ~ ">5 years"
  )) %>%
  mutate(subjectResponse=annotateSubjectResponse(Subject)) %>%
  ggplot() +
  geom_bar(aes(x=fct_relevel(gsub(" ","\n",LastAbxDisplay),
                             gsub(" ", "\n", c("<6 months", "6-12 months", "1-2 years",
                               "3-5 years", ">5 years", "Unknown"))),
               fill=fct_relevel(subjectResponse, rev(ABXRESPONSE)))) +
  scale_fill_manual(values=PALETTE.ABXRESPONSE, name="Subject\nresponse") +
  xlab("Most recent antibiotic usage") + ylab("Number of\nparticipants") +
  guides(fill="none") +
  DEFAULTS.THEME_PRINT + theme(axis.text.x=element_text(size=5))
savePNGPDF(paste0(OUTDIR, "3-subjectResponse-priorAbx"), p3_recentAbx, suppfig3height, 2.25)

# Calculate the median diversity in numbers of effective species before antibiotics.
dataPreAbxDiversity <- dataAlphaEffectiveSpecies %>%
  group_by(subject, subjectResponse) %>%
  summarize(medianPreAbxDiversity=median(value[timepoint<=29]))
p3_preAbxDiversity <- dataPreAbxDiversity %>%
  filter(subjectResponse!="control") %>%
  ggplot() +
  geom_point(aes(x=fct_relevel(gsub(" ","\n",subjectResponse), names(PALETTE.ABXRESPONSELINEBREAK)), 
                 y=medianPreAbxDiversity,
                 color=fct_relevel(subjectResponse, ABXRESPONSE)),
             position=position_jitter(width=0.2, height=0, seed=0), size=1, alpha=0.8) +
  scale_color_manual(values=PALETTE.ABXRESPONSE, name="Subject\nresponse") + ylim(0,NA) +
  xlab("Subject response") + ylab("Pre-antibiotic\ndiversity") +
  DEFAULTS.THEME_PRINT + theme(axis.text.x=element_text(size=5)) +
  guides(color = guide_legend(override.aes = list(size = 0.75))) +
  theme(legend.key.size = unit(0.5, "lines"))
savePNGPDF(paste0(OUTDIR, "3-subjectResponse-preAbxDiversity"), 
           p3_preAbxDiversity + guides(color="none"),
           suppfig3height, 1.25)
# No significant relationship between pre-abx diversity and subject response (categorical)
# Wilcoxon p=0.22
wilcox.test(dataPreAbxDiversity %>% filter(subjectResponse %in% c("transient response")) %>% pull(medianPreAbxDiversity),
            dataPreAbxDiversity %>% filter(subjectResponse %in% c("lasting response")) %>% pull(medianPreAbxDiversity))

# Import the CFU data.
dataCFUs <- read.table("workflow/analysis/calculateSampleCFUs/out/e0025-CFUcalculations.tsv",
                       header=TRUE, stringsAsFactors = FALSE)
# Parse sample metadata.
dataCFUs <- dataCFUs %>%
  mutate(subject=substr(sample,1,3), timepoint=as.numeric(substr(sample,5,7)))
# Summarize CFU changes.
deltaCFUs <- dataCFUs %>%
  filter(dilutionReplicate==1, plateReplicate==1) %>%
  group_by(subject) %>% arrange(timepoint) %>%
  mutate(timePeriodIndex=row_number()) %>%
  mutate(timePeriod=case_when(
    timePeriodIndex==1 ~ "initial",
    timePeriodIndex==2 ~ "preAbx",
    timePeriodIndex==3 ~ "postAbx",
    timePeriodIndex==4 ~ "final"
  )) %>% ungroup() %>%
  dplyr::select(subject, timePeriod, CFUsPerGStool) %>%
  pivot_wider(names_from=timePeriod, values_from=CFUsPerGStool) %>%
  mutate(initialPreAbxDeltaCFU=log10(preAbx)-log10(initial),
         prePostAbxDeltaCFU=log10(postAbx)-log10(preAbx),
         postAbxFinalDeltaCFU=log10(final)-log10(postAbx),
         initialPostAbxDeltaCFU=log10(postAbx)-log10(initial),
         initialFinalDeltaCFU=log10(final)-log10(initial)) %>%
  mutate(subjectResponse=annotateSubjectResponse(subject)) %>%
  filter(subjectResponse!="control") %>%
  dplyr::select(subject, subjectResponse, contains("DeltaCFU")) %>%
  pivot_longer(contains("DeltaCFU"), names_to="timePeriod", values_to="deltaCFU") %>%
  mutate(timePeriodDisplay= case_when(
    timePeriod=="initialPreAbxDeltaCFU" ~ "pre-abx\nday 0",
    timePeriod=="prePostAbxDeltaCFU" ~ "days 0 to 7\nabx",
    timePeriod=="postAbxFinalDeltaCFU" ~ "days 7 to 35\npost-abx",
    timePeriod=="initialPostAbxDeltaCFU" ~ "post-abx\nday 7",
    timePeriod=="initialFinalDeltaCFU" ~ "final\nday 35"
  ))
p3_deltaCFUs <- deltaCFUs %>%
  filter(timePeriod %in% c("initialPreAbxDeltaCFU","initialPostAbxDeltaCFU","initialFinalDeltaCFU")) %>%
  arrange(desc(subjectResponse)) %>%
  ggplot() +
  geom_point(aes(x=fct_relevel(timePeriodDisplay, 
                               c("pre-abx\nday 0", "post-abx\nday 7", "final\nday 35")), 
                 y=deltaCFU,
                 color=fct_relevel(subjectResponse, ABXRESPONSE)), alpha=0.8, size=1,
             position=position_jitter(width=0.2, height=0, seed=0)) +
  scale_color_manual(values=PALETTE.ABXRESPONSE) + guides(color="none") +
  xlab("Time period") + ylab("log10 change\nfrom study start,\nCFUs per g stool") +
  DEFAULTS.THEME_PRINT + theme(axis.title.x=element_blank())
savePNGPDF(paste0(OUTDIR, "3-deltaCFUs"), p3_deltaCFUs, suppfig3height, 2)


# Import the proportion of the gut microbiome that consisted of shared strains
# at the beginning of the study.
# Plot the relative abundance of shared species and strains
# at the initial timepoint for subjects in study arm X.
# Import the relative abundance data annotated by strain sharing.
# Plot the relative abundance of shared species and strains
# at the initial timepoint for subjects in study arm X.
# Import the summary of the relative abundance data.
dataRelAbundanceStrainSharingSummaryRaw <-
  read.table(gzfile("workflow/analysis/integrateSpeciesAbundancesSharedStrains/out/annotatedSpeciesAbundances-majorTimepoints-summary.txt"),
             header=TRUE, stringsAsFactors = FALSE)
# Also use complete to infer when certain subjects pairs have no strains in common.
dataRelAbundanceStrainSharingSummary <- dataRelAbundanceStrainSharingSummaryRaw %>%
  filter(sample %in% samplesInitial, annotationSample %in% samplesInitial,
         sample!=annotationSample, sample %in% samplesX, annotationSample %in% samplesX,
         substr(sample,1,3)!=substr(annotationSample,1,3)) %>%
  mutate(samplePair=paste0(sample,"_", annotationSample)) %>%
  complete(samplePair, annotation, fill = list(totRelAbundance=0, numSpecies=0)) %>%
  mutate(sample=substr(samplePair,1,7), annotationSample=substr(samplePair,9,15)) %>%
  filter(annotation=="strainsShared") %>%
  group_by(annotation, sample) %>% summarize(maxRelAbundanceShared=max(totRelAbundance)) %>%
  mutate(subject=substr(sample,1,3)) %>%
  mutate(subjectResponse=annotateSubjectResponse(subject))
# Summarize the proportion of the gut microbiome that consists of shared strains
# between cohabiting subjects.
p3_strainSharing <- dataRelAbundanceStrainSharingSummary %>%
  filter(subjectResponse!="control") %>%
  ggplot() +
  geom_point(aes(x=fct_relevel(gsub(" ","\n",subjectResponse), names(PALETTE.ABXRESPONSELINEBREAK)), 
                 y=maxRelAbundanceShared,
                 color=factor(subjectResponse)),
             position=position_jitter(width=0.2, height=0, seed=0), size=1, alpha=0.8) +
  scale_color_manual(values=PALETTE.ABXRESPONSE, name="Subject\nresponse") +
  xlab("Subject response") + ylab("Total abundance,\nshared strains") +
  ylim(0,1) + guides(color="none") +
  DEFAULTS.THEME_PRINT + theme(axis.text.x=element_text(size=5))
savePNGPDF(paste0(OUTDIR,"3-strainSharing-abxResponses"), p3_strainSharing,
           suppfig3height, suppfig3height)
# Statistically compare the proportion of the microbiome that consists of shared strains
# in subjects with transient and lasting responses.
# Wilcoxon rank-sum p=0.8753 (NOTE that there are ties)
wilcox.test(dataRelAbundanceStrainSharingSummary %>% filter(subjectResponse=="transient response") %>% pull(maxRelAbundanceShared),
            dataRelAbundanceStrainSharingSummary %>% filter(subjectResponse=="lasting response") %>% pull(maxRelAbundanceShared))

# Import species abundances and defaults for relative abundance plots.
source("workflow/analysis/generateSpeciesAbundances/loadSpeciesAbundancePlots.R")
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
subjectsOrderedJSD <- dataBetaJSDMaxFinal %>%
  arrange(maxJSD) %>% pull(subject)

p3_allRelativeAbundancePlots <- plotSpeciesAbundanceDaysFromAbxStart(speciesAbundances %>%
                                       filter(sample %in% samplesXmain, subject %in% subjectsAbx) %>%
                                       group_by(subject, species_id) %>% filter(max(relative_abundance)>LIMITOFDETECTION)) +
  facet_wrap(~fct_relevel(subject, subjectsOrderedJSD), scales="free") +
  DEFAULTS.THEME_PRINT
savePNGPDF(paste0(OUTDIR, "3-relativeAbundancesAllAbx"), p3_allRelativeAbundancePlots,
           5, 6.8)

# Figures 4-7 - heatmaps, subjects with lasting responses -----------------


source("workflow/analysis/plotHelpers-heatmaps.R")
OUTDIR <- "workflow/analysis/supp-figures-v2/out/"

# Import species abundances and defaults for relative abundance plots.
source("workflow/analysis/generateSpeciesAbundances/loadSpeciesAbundancePlots.R")

plotSpeciesAbundanceByTimepoint <- function(x){
  x %>%
    ungroup() %>% group_by(subject) %>%
    mutate(timepoint=timepoint-29, minTimepoint=min(timepoint), maxTimepoint=max(timepoint)) %>%
    ggplot() +
    geom_area(aes(x=factor(timepoint), y=relative_abundance,
                  group=fct_relevel(species_id, commonSpecies),
                  fill=factor(color), alpha=alpha),
              color="black", linewidth=0.1) +
    facet_grid(.~(subject), scales="free_x", space="free_x") +
    scale_fill_manual(values=commonSpeciesPalette) +
    scale_alpha_identity() +
    ylim(0,1) +
    guides(fill="none") +
    xlab("Study day") + ylab("Relative\nabundance") +
    DEFAULTS.THEME_PRINT + 
    theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=4.5),
          axis.text.y=element_text(size=4.5), strip.text.x=element_text(size=5))
}

speciesAbundancesFiltered <- speciesAbundances %>%
  group_by(species_id) %>% filter(max(relative_abundance)>LIMITOFDETECTION)

# Subject XAA: arrange the heatmap and relative abundance plots on the same grid.
p4_arranged <- 
  (plotSpeciesAbundanceByTimepoint(speciesAbundancesFiltered %>% filter(subject %in% c("XAA","XAC","XAD"))) + 
     theme(axis.title.y=element_blank())) / (plotHeatmap("XAA") + guides(fill="none", color="none", shape="none")) +
  plot_layout(heights=c(1, 7))
p4_arranged
savePNGPDF(paste0(OUTDIR, "4-relAbundance-heatmap-XAA"), p4_arranged, 7, 6)
savePNGPDF(paste0(OUTDIR, "4-heatmap-legend-XAA"), get_legend(plotHeatmap("XAA")), 2.5, 0.75)

# Subject XBA: arrange the heatmap and relative abundance plots on the same grid.
p5_arranged <- 
  (plotSpeciesAbundanceByTimepoint(speciesAbundancesFiltered %>% 
                                     filter(subject %in% c("XBA","XBB"), 
                                            !(sample %in% c(samplesXBAextraMain, samplesXBAextraFollowup)))) + 
     theme(axis.title.y=element_blank())) / (plotHeatmap("XBA") + guides(fill="none", color="none", shape="none")) +
  plot_layout(heights=c(1, 6))
p5_arranged
savePNGPDF(paste0(OUTDIR, "5-relAbundance-heatmap-XBA"), p5_arranged, 7, 6.5)
savePNGPDF(paste0(OUTDIR, "5-heatmap-legend-XBA"), get_legend(plotHeatmap("XBA")), 2.5, 0.75)

# Subject XEA: arrange the heatmap and relative abundance plots on the same grid.
p6_arranged <- 
  (plotSpeciesAbundanceByTimepoint(speciesAbundancesFiltered %>% 
                                     filter(subject %in% c("XEA","XEB"))) + 
     theme(axis.title.y=element_blank())) / (plotHeatmap("XEA") + guides(fill="none", color="none", shape="none")) +
  plot_layout(heights=c(1, 6))
p6_arranged
savePNGPDF(paste0(OUTDIR, "6-relAbundance-heatmap-XEA"), p6_arranged, 6, 4)
savePNGPDF(paste0(OUTDIR, "6-heatmap-legend-XEA"), get_legend(plotHeatmap("XEA")), 2.5, 0.75)

# Subject XKA: arrange the heatmap and relative abundance plots on the same grid.
p7_arranged <- 
  (plotSpeciesAbundanceByTimepoint(speciesAbundancesFiltered %>% 
                                     filter(subject %in% c("XKA","XKB"))) + 
     theme(axis.title.y=element_blank())) / (plotHeatmap("XKA") + guides(fill="none", color="none", shape="none")) +
  plot_layout(heights=c(1, 5.6))
p7_arranged
savePNGPDF(paste0(OUTDIR, "7-relAbundance-heatmap-XKA"), p7_arranged, 5, 5)
savePNGPDF(paste0(OUTDIR, "7-heatmap-legend-XKA"), get_legend(plotHeatmap("XKA")), 2.5, 0.75)

# Generate the legend for the species coloring.
PALETTE.SPECIESTRAJECTORIES
speciesTrajectoryLegend <- ggplot() +
  geom_tile(aes(x=names(PALETTE.SPECIESTRAJECTORIES)[-2], y=seq(1:4),
                fill=names(PALETTE.SPECIESTRAJECTORIES)[-2])) +
  scale_fill_manual(values=PALETTE.SPECIESTRAJECTORIES[-2], name="Species\ntrajectories",
                    breaks=c("notDisrupted","recovered","notRecovered","colonized"),
                    labels=c("species\nmaintained","species\ndisrupted, recovered",
                             "species\ndisrupted, not recovered","species\ncolonized")) +
  DEFAULTS.THEME_PRINT +
  guides(fill = guide_legend(override.aes = list(size = 0.5))) +
  theme(legend.key.size = unit(0.5, "lines"))
savePNGPDF(paste0(OUTDIR, "4-speciesLegend"), get_legend(speciesTrajectoryLegend), 1, 1.5)

# Export the common families palette.
# List only the top 8 families.
topFamilies <- (speciesAbundances %>%
                  filter(sample %in% samplesInitial) %>%
                  group_by(family) %>% summarize(totAbundance=sum(relative_abundance)) %>%
                  top_n(10, totAbundance) %>%
                  ungroup() %>% mutate(family=ifelse(family=="","Unknown",family)))$family
topFamiliesPalette <- commonFamiliesPalette[topFamilies]
topFamiliesdataframe <- data.frame(family=topFamilies, color=topFamiliesPalette)
topFamiliesdataframe$family <- factor(topFamilies,
                                      levels=c("Unknown", topFamilies[-which(topFamilies=="Unknown")]))
p4_familiesPalette <- topFamiliesdataframe %>%
  ggplot() +
  geom_tile(aes(x=family, y=0, fill=factor(family))) +
  scale_fill_manual(name="Family", values=topFamiliesPalette) +
  DEFAULTS.THEME_PRINT +
  guides(fill = guide_legend(override.aes = list(size = 0.5))) +
  theme(legend.key.size = unit(0.5, "lines"))
savePNGPDF(paste0(OUTDIR, "4-relativeAbundance-palette"), get_legend(p4_familiesPalette),
           1.5, 1)
savePNGPDF(paste0(OUTDIR, "4-relativeAbundance-palette-small"), 
           get_legend(p4_familiesPalette + theme(legend.text=element_text(size=4.5)) +
             guides(fill = guide_legend(override.aes = list(size = 0.3))) +
             theme(legend.key.size = unit(0.3, "lines"))),
           1.5, 1)

# Figure 8 - species and strain dynamics, main study ----------------------

source("workflow/analysis/plotHelpers.R")
OUTDIR <- "workflow/analysis/supp-figures-v2/out/"

suppfig8height <- 1.25

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
subjectsOrderedJSD <- dataBetaJSDMaxFinal %>%
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
  )) %>% mutate(speciesTrajectory=fct_relevel(speciesTrajectory, c("colonized", "notRecovered","recovered","notDisrupted")))
p8_numSpeciesPerCategory_controls <- dataNumSpeciesPerCategory %>%
  mutate(subjectResponse=annotateSubjectResponse(subject),
         subjectResponseDisplay=annotateSubjectResponseNumSubjects(subject)) %>%
  filter(subjectResponse=="control") %>%
  ggplot() +
  geom_bar(aes(x=fct_relevel(subject, subjectsOrderedJSD), y=numSpecies, 
               fill=fct_relevel(speciesTrajectory, rev(names(PALETTE.SPECIESTRAJECTORIES)))), 
           stat="identity") +
  facet_grid(.~fct_rev(gsub(" \\(", "\n\\(", subjectResponseDisplay)), scales="free_x", space="free_x") +
  scale_fill_manual(values=PALETTE.SPECIESTRAJECTORIES, name="Species trajectory",
                    labels=c("new species","species disrupted, not recovered",
                             "species disrupted, recovered", "species maintained")) +
  xlab("Subject") + ylab("Number of\nspecies") +
  guides(fill = guide_legend(override.aes = list(size = 1))) +
  theme(legend.key.size = unit(0.5, "lines")) +
  DEFAULTS.THEME_PRINT + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=5))
p8_numSpeciesPerCategory_controls
savePNGPDF(paste0(OUTDIR,"8-numSpeciesPerSubject-perAnnotation"),
           p8_numSpeciesPerCategory_controls, 1.5, 5.25)

# Calculate the total proportion of the microbiome at the initial and final timepoints
# that consists of disrupted, recovered, and colonized species.
dataSpeciesTrajectories <- read.table("workflow/analysis/annotateSpeciesTrajectoriesSharedStrains/out/speciesAbundances-trajectoriesSharingAnnotated-all.txt.gz",
                                      header=TRUE, stringsAsFactors = FALSE)
dataSpeciesTrajectoriesRelAbundance <- dataSpeciesTrajectories %>%
  mutate(speciesTrajectory=case_when(
    speciesStatusAbxAllTimepoints=="speciesNotDisrupted" ~ "notDisrupted",
    speciesStatusAbxAllTimepoints=="speciesRecoveredMain" ~ "recovered",
    speciesStatusAbxAllTimepoints=="speciesRecoveredFollowup" | speciesStatusAbxAllTimepoints=="speciesNotRecovered" ~ "notRecovered",
    speciesStatusColonizationAllTimepoints!="speciesNotColonized" | (strainTurnover & timepoint>timeOfStrainTurnover) ~ "colonized"
  )) %>%
  filter(sample %in% c(samplesInitial, samplesFinal)) %>%
  mutate(timePeriod=ifelse(sample %in% samplesInitial, "initial", "final"),
         speciesTrajectory=ifelse(is.na(speciesTrajectory), "notDisrupted", speciesTrajectory)) %>%
  group_by(subject, timePeriod, speciesTrajectory) %>%
  summarize(totAbundance=sum(relative_abundance)) %>%
  mutate(subjectResponse=annotateSubjectResponse(subject)) 
dataSpeciesTrajectoriesRelAbundance <- dataSpeciesTrajectoriesRelAbundance %>%
  pivot_wider(names_from="speciesTrajectory", values_from="totAbundance", values_fill=0) %>%
  mutate(notClassified=1-colonized-notDisrupted-notRecovered-recovered) %>%
  pivot_longer(colonized:notClassified, names_to="speciesTrajectory", values_to="totAbundance")
p8_speciesTrajectory_totAbundance <- dataSpeciesTrajectoriesRelAbundance %>%
  mutate(subjectResponseDisplay=annotateSubjectResponseNumSubjects(subject)) %>%
  mutate(timePeriodDisplay=ifelse(timePeriod=="initial","day -28","day 35")) %>%
  ggplot() +
  geom_bar(aes(x=fct_relevel(subject, subjectsOrderedJSD), y=totAbundance,
               fill=fct_relevel(speciesTrajectory, rev(names(PALETTE.SPECIESTRAJECTORIES)))), 
           stat="identity") +
  facet_grid((timePeriodDisplay)~
               fct_relevel(gsub(" \\(", "\n\\(", subjectResponseDisplay), gsub(" \\(", "\n\\(",subjectResponseNumSubjectsOrdered)), 
             scales="free", space="free") +
  scale_fill_manual(values=PALETTE.SPECIESTRAJECTORIES, name="Species trajectory",
                    labels=c("new species","species disrupted, not recovered",
                             "species disrupted, recovered", "species not classified", "species maintained")) +
  xlab("Subject") + ylab("Total abundance") +
  guides(fill = guide_legend(override.aes = list(size = 1))) +
  theme(legend.key.size = unit(0.5, "lines")) +
  DEFAULTS.THEME_PRINT + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=5)) +
  theme(strip.text.x=element_text(size=6))
savePNGPDF(paste0(OUTDIR, "8-totalAbundancePerSubject-annotation"),
           p8_speciesTrajectory_totAbundance + guides(fill="none"), 2, 6.9)  


# Plot the taxonomic enrichment of species that are disrupted and do not recover.
# Import the data on taxonomic enrichment among disrupted species that do not recover.
dataEnrichmentDisruptedNoRecovery <- 
  read.table("workflow/analysis/calculateTaxonomicEnrichment/out/disruptedSpeciesNoRecovery-randomDistributionSummary.txt",
             header=TRUE, stringsAsFactors = FALSE, sep="\t")
# Plot a volcano plot of the family enrichment analysis.
dataEnrichmentDisruptedNoRecovery <- dataEnrichmentDisruptedNoRecovery %>%
  mutate(log2Enrichment=log2((meanActual+1)/(meanRandom+1)),
         significant=(pvalue<lowerThreshold),
         pvalue=ifelse(pvalue==0, 1/numCounts, pvalue)) 
p8_enrichment_disrupted_noRecovery <- dataEnrichmentDisruptedNoRecovery %>%
  ggplot() +
  geom_hline(aes(yintercept=-log10(lowerThreshold)), linetype="dashed", alpha=0.5) +
  geom_point(aes(x=log2Enrichment, y=-log10(pvalue),
                 color=factor(ifelse(!significant, "ns",
                                     ifelse(log2Enrichment>0,"sig-plus","sig-minus")))), 
             alpha=0.8, size=1) +
  scale_color_manual(values=c("gray40","dodgerblue3","firebrick")) +
  guides(color="none") + ylim(0,4.5) + xlim(-2.5,2.5) +
  xlab("Fold enrichment,\nno recovery") + 
  ylab("Significance\n-log10(p-value)") +
  DEFAULTS.THEME_PRINT
savePNGPDF(paste0(OUTDIR, "8-enrichment-volcano-disruptedNoRecovery"), 
           p8_enrichment_disrupted_noRecovery, 
           suppfig8height, 1.5)

# Calculate the proportion of each subject's gut microbiome at the initial timepoint
# that is comprised of taxa that are more or less susceptible to disruption.
# Import the data on taxonomic enrichment among disrupted species.
dataEnrichmentDisrupted <- 
  read.table("workflow/analysis/calculateTaxonomicEnrichment/out/disruptedSpecies-randomDistributionSummary.txt",
             header=TRUE, stringsAsFactors = FALSE, sep="\t")
# Calculate the direction and amount of enrichment.
dataEnrichmentDisrupted <- dataEnrichmentDisrupted %>%
  mutate(log2Enrichment=log2((meanActual+1)/(meanRandom+1)),
         significant=(pvalue<lowerThreshold),
         pvalue=ifelse(pvalue==0, 1/numCounts, pvalue)) 
taxaEnrichedDisrupted <- dataEnrichmentDisrupted %>% filter(log2Enrichment>0, significant) %>% pull(group)
taxaDepletedDisrupted <- dataEnrichmentDisrupted %>% filter(log2Enrichment<0, significant) %>% pull(group)
# Import each subject's gut microbiome composition at the initial timepoint.
source("workflow/analysis/generateSpeciesAbundances/loadSpeciesAbundancesFiltered.R")
# Analyze only community composition at the initial timepoint.
dataInitialCompositionX <- speciesAbundances %>%
  filter(sample %in% samplesInitial, subject %in% subjectsX)
dataInitialCompositionXTaxaDisruption <- dataInitialCompositionX %>%
  mutate(annotation=case_when(
    family %in% taxaEnrichedDisrupted ~ "more likely",
    family %in% taxaDepletedDisrupted ~ "less likely",
    !(family %in% c(taxaEnrichedDisrupted, taxaDepletedDisrupted)) ~ "none"
  )) %>%
  group_by(subject, annotation) %>% summarize(totRelAbundance=sum(relative_abundance)) %>%
  mutate(subjectResponse=annotateSubjectResponse(subject))
p8_abundanceTaxonomicEnrichmentDisruption <- dataInitialCompositionXTaxaDisruption %>%
  filter(subject %in% subjectsAbx, annotation!="none") %>%
  ggplot() +
  geom_point(aes(x=fct_relevel(gsub(" ","\n",subjectResponse), names(PALETTE.ABXRESPONSELINEBREAK)), 
                 y=totRelAbundance, color=subjectResponse), 
             alpha=0.8, size=1, position=position_jitter(width=0.2, height=0, seed=0)) +
  facet_wrap(~annotation) + ylim(0,NA) +
  scale_color_manual(values=PALETTE.ABXRESPONSE) + guides(color="none") +
  ylab("Total abundance") + xlab("Subject response") + ylim(0,1) +
  DEFAULTS.THEME_PRINT + theme(axis.text.x=element_text(size=5))
savePNGPDF(paste0(OUTDIR, "8-taxEnrichmentAbundance-disruption"), 
           p8_abundanceTaxonomicEnrichmentDisruption, suppfig8height, 2.5)
# Statistically compare the proportion of more and less easily disrupted taxa.
# total abundance of taxa that are more likely to be disrupted
# Wilcoxon rank-sum p=0.14 (not significant)
wilcox.test(dataInitialCompositionXTaxaDisruption %>%
              filter(subjectResponse=="transient response", annotation=="more likely") %>% pull(totRelAbundance),
            dataInitialCompositionXTaxaDisruption %>%
              filter(subjectResponse=="lasting response", annotation=="more likely") %>% pull(totRelAbundance))
# Statistically compare the proportion of more and less easily disrupted taxa.
# total abundance of taxa that are less likely to be disrupted
# Wilcoxon rank-sum p=0.19 (not significant)
wilcox.test(dataInitialCompositionXTaxaDisruption %>%
              filter(subjectResponse=="transient response", annotation=="less likely") %>% pull(totRelAbundance),
            dataInitialCompositionXTaxaDisruption %>%
              filter(subjectResponse=="lasting response", annotation=="less likely") %>% pull(totRelAbundance))


# Calculate the proportion of each subject's gut microbiome at the initial timepoint
# that is comprised of taxa that are more or less susceptible to disruption without recovery.
# Import the data on taxonomic enrichment among disrupted species that do not recover.
dataEnrichmentDisruptedNoRecovery <- 
  read.table("workflow/analysis/calculateTaxonomicEnrichment/out/disruptedSpeciesNoRecovery-randomDistributionSummary.txt",
             header=TRUE, stringsAsFactors = FALSE, sep="\t")
# Plot a volcano plot of the family enrichment analysis.
dataEnrichmentDisruptedNoRecovery <- dataEnrichmentDisruptedNoRecovery %>%
  mutate(log2Enrichment=log2((meanActual+1)/(meanRandom+1)),
         significant=(pvalue<lowerThreshold),
         pvalue=ifelse(pvalue==0, 1/numCounts, pvalue)) 
taxaEnrichedDisruptedNoRecovery <- dataEnrichmentDisruptedNoRecovery %>% filter(log2Enrichment>0, significant) %>% pull(group)
taxaDepletedDisruptedNoRecovery <- dataEnrichmentDisruptedNoRecovery %>% filter(log2Enrichment<0, significant) %>% pull(group)
# Import each subject's gut microbiome composition at the initial timepoint.
source("workflow/analysis/generateSpeciesAbundances/loadSpeciesAbundancesFiltered.R")
# Analyze only community composition at the initial timepoint.
dataInitialCompositionX <- speciesAbundances %>%
  filter(sample %in% samplesInitial, subject %in% subjectsX)
dataInitialCompositionXTaxaDisruptionNoRecovery <- dataInitialCompositionX %>%
  mutate(annotation=case_when(
    family %in% taxaEnrichedDisruptedNoRecovery ~ "more likely",
    family %in% taxaDepletedDisruptedNoRecovery ~ "less likely",
    !(family %in% c(taxaEnrichedDisruptedNoRecovery, taxaDepletedDisruptedNoRecovery)) ~ "none"
  )) %>%
  group_by(subject, annotation) %>% summarize(totRelAbundance=sum(relative_abundance)) %>%
  mutate(subjectResponse=annotateSubjectResponse(subject))
p8_abundanceTaxonomicEnrichmentDisruptionNoRecovery <- dataInitialCompositionXTaxaDisruptionNoRecovery %>%
  filter(subject %in% subjectsAbx, annotation!="none") %>%
  ggplot() +
  geom_point(aes(x=fct_relevel(gsub(" ","\n",subjectResponse), names(PALETTE.ABXRESPONSELINEBREAK)), 
                 y=totRelAbundance, color=subjectResponse), 
             alpha=0.8, size=1, position=position_jitter(width=0.2, height=0, seed=0)) +
  facet_wrap(~annotation) + ylim(0,NA) +
  scale_color_manual(values=PALETTE.ABXRESPONSE) + guides(color="none") +
  ylab("Total abundance") + xlab("Subject response") + ylim(0,1) +
  DEFAULTS.THEME_PRINT + theme(axis.text.x=element_text(size=5))
savePNGPDF(paste0(OUTDIR, "8-taxEnrichmentAbundance-disruptionNoRecovery"), 
           p8_abundanceTaxonomicEnrichmentDisruptionNoRecovery, suppfig8height, 2.5)
# Statistically compare the proportion of more and less easily disrupted taxa.
# total abundance of taxa that are more likely to be disrupted
# Wilcoxon rank-sum p=0.67 (not significant)
wilcox.test(dataInitialCompositionXTaxaDisruptionNoRecovery %>%
              filter(subjectResponse=="transient response", annotation=="more likely") %>% pull(totRelAbundance),
            dataInitialCompositionXTaxaDisruptionNoRecovery %>%
              filter(subjectResponse=="lasting response", annotation=="more likely") %>% pull(totRelAbundance))
# Statistically compare the proportion of more and less easily disrupted taxa.
# total abundance of taxa that are less likely to be disrupted
# Wilcoxon rank-sum p=0.22 (not significant)
wilcox.test(dataInitialCompositionXTaxaDisruptionNoRecovery %>%
              filter(subjectResponse=="transient response", annotation=="less likely") %>% pull(totRelAbundance),
            dataInitialCompositionXTaxaDisruptionNoRecovery %>%
              filter(subjectResponse=="lasting response", annotation=="less likely") %>% pull(totRelAbundance))

# Import species trajectory summary data.
dataSpeciesTrajectoriesSummary <- 
  read.table("workflow/analysis/annotateSpeciesTrajectoriesSharedStrains/out/speciesTrajectoriesSharingSummary.txt",
             header=TRUE, stringsAsFactors = FALSE)
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
# strains shared, species recovered: 55
# strains shared, species not recovered: 14
# strains not shared, species recovered: 106
# strains not shared, species not recovered: 35
dataSharingRecoveryStrainsMinimalTransient <- dataSharingRecoveryMinimalTransient %>%
  group_by(speciesRecoveredMain, strainSharing) %>%
  summarize(numSpecies=n()) %>%
  filter(!is.na(strainSharing)) %>%
  group_by(strainSharing) %>%
  mutate(pctSpecies=numSpecies/sum(numSpecies),
         totalSpecies=sum(numSpecies))
# chi-square p=0.5784
sharingRecoveryStrainsMinimalTransient_chisq <- 
  chisq.test(matrix(dataSharingRecoveryStrainsMinimalTransient$numSpecies, nrow=2, ncol=2))


# Import the annotated species trajectories to plot the number of new colonizers over time.
dataSpeciesTrajectories <-
  read.table("workflow/analysis/annotateSpeciesTrajectoriesSharedStrains/out/speciesAbundances-trajectoriesSharingAnnotated-all.txt.gz",
             header=TRUE, stringsAsFactors = FALSE)

# Plot the number of colonizing strains in each antibiotic-taking subject.
dataColonizingSpeciesMain <- dataSpeciesTrajectories %>%
  mutate(speciesTrajectory=ifelse((timeOfColonizationAllTimepoints<75 | timeOfStrainTurnover<75),
                                  "colonizedStrainTurnover", speciesTrajectory)) %>%
  group_by(subject, speciesTrajectory) %>% 
  summarize(numStrains=n_distinct(species_id)) %>% ungroup() %>%
  complete(subject, speciesTrajectory, fill=list(numStrains=0)) %>%
  filter(speciesTrajectory=="colonizedStrainTurnover") %>%
  mutate(subjectTreatment=annotateSubjectTreatment(subject),
         subjectResponse=annotateSubjectResponse(subject),
         subjectResponseDisplay=gsub(" ","\n",subjectResponse)) %>%
  arrange(desc(numStrains))
p8_numColonizingStrains <- dataColonizingSpeciesMain %>%
  ggplot() +
  geom_point(aes(x=factor(subjectResponseDisplay, names(PALETTE.ABXRESPONSELINEBREAK)), 
                 y=numStrains, color=factor(subjectResponseDisplay)),
             position=position_jitter(width=0.2, height=0, seed=0), alpha=0.8, size=1) +
  scale_color_manual(values=PALETTE.ABXRESPONSELINEBREAK) + guides(color="none") +
  xlab("Subject response") + ylab("Number of colonizers,\nfirst month post-abx") +
  DEFAULTS.THEME_PRINT + theme(axis.title.y=element_text(size=6))
savePNGPDF(paste0(OUTDIR, "8-numColonizingStrains-main"), p8_numColonizingStrains,
           suppfig8height, 1.9)
# Statistically compare the number of colonizers for each subject response group.
# control vs. minimal/transient: p=0.01839 (Wilcoxon)
# control vs. lasting: p=3.892e-5
# minimal/transient vs. lasting: p=0.001795
wilcox.test(dataColonizingSpeciesMain %>% filter(subjectResponse=="control") %>% pull(numStrains),
            dataColonizingSpeciesMain %>% filter(subjectResponse=="transient response") %>% pull(numStrains))
wilcox.test(dataColonizingSpeciesMain %>% filter(subjectResponse=="control") %>% pull(numStrains),
            dataColonizingSpeciesMain %>% filter(subjectResponse=="lasting response") %>% pull(numStrains))
wilcox.test(dataColonizingSpeciesMain %>% filter(subjectResponse=="transient response") %>% pull(numStrains),
            dataColonizingSpeciesMain %>% filter(subjectResponse=="lasting response") %>% pull(numStrains))


# Import the total abundance of species over time,
# annotated based on all timepoints.
dataTotalDisruptedColonized_allTimepoints <-
  read.table("workflow/analysis/annotateSpeciesTrajectoriesSharedStrains/out/speciesTrajectoriesPerSubject-totalAbundance-allTimepoints.txt",
             header=TRUE, stringsAsFactors = FALSE)
dataTotalDisruptedColonized_allTimepointsSummary <- dataTotalDisruptedColonized_allTimepoints %>%
  filter(sample %in% samplesFinal, speciesTrajectory %in% c("colonized", "strainTurnover")) %>%
  group_by(subject) %>% 
  summarize(speciesTrajectory="colonizedStrainTurnover", totalAbundance=sum(totalAbundance)) %>%
  mutate(subjectTreatment=annotateSubjectTreatment(subject),
         subjectResponse=annotateSubjectResponse(subject),
         subjectResponseDisplay=gsub(" ","\n",subjectResponse))
# Plot the proportion of the gut microbiome annotated as newly colonizing over time.
p8_pctColonizedAll <- dataTotalDisruptedColonized_allTimepointsSummary %>%
  ggplot() +
  geom_point(aes(x=factor(subjectResponseDisplay, names(PALETTE.ABXRESPONSELINEBREAK)), 
                 y=totalAbundance, color=factor(subjectResponseDisplay)),
             position=position_jitter(width=0.2, height=0, seed=0), alpha=0.8, size=1) +
  scale_color_manual(values=PALETTE.ABXRESPONSELINEBREAK) + guides(color="none") +
  xlab("Subject response") + ylab("Total abundance of\nnew colonizers,\nfirst month post-abx") +
  DEFAULTS.THEME_PRINT + theme(axis.title.y=element_text(size=6)) + ylim(0,1)
p8_pctColonizedAll
savePNGPDF(paste0(OUTDIR, "8-pctColonizedTurnover-byResponse"), p8_pctColonizedAll, 
           suppfig8height, 1.9)
# Statistically assess the proportion of the gut microbiome comprised of new colonizers
# at the end of the main study.
# Wilcoxon, control vs. minimal/transient: p=0.09871
# Wilcoxon, control vs. lasting: p=0.0001873
# Wilcoxon, minimal/transient vs. lasting: p=0.01189
wilcox.test(dataTotalDisruptedColonized_allTimepointsSummary %>% filter(subjectResponse=="control") %>% pull(totalAbundance),
            dataTotalDisruptedColonized_allTimepointsSummary %>% filter(subjectResponse=="transient response") %>% pull(totalAbundance))
wilcox.test(dataTotalDisruptedColonized_allTimepointsSummary %>% filter(subjectResponse=="control") %>% pull(totalAbundance),
            dataTotalDisruptedColonized_allTimepointsSummary %>% filter(subjectResponse=="lasting response") %>% pull(totalAbundance))
wilcox.test(dataTotalDisruptedColonized_allTimepointsSummary %>% filter(subjectResponse=="transient response") %>% pull(totalAbundance),
            dataTotalDisruptedColonized_allTimepointsSummary %>% filter(subjectResponse=="lasting response") %>% pull(totalAbundance))


# Import the HMP prevalence and abundance of resident species and new colonizers.
residentSpeciesColonizersHMPprevalenceAbundance <- 
  read.table("workflow/analysis/calculateHMPspeciesPrevalenceAbundance/out/residentColonizingSpecies-prevalenceAbundance.txt",
             header=TRUE, stringsAsFactors = FALSE) %>%
  mutate(type=fct_relevel(type, c("resident\nspecies", "new colonizer\n<1 month post-abx", "new colonizer\n>1 month post-abx")))
p8_medianAbundanceResidentColonizers <- residentSpeciesColonizersHMPprevalenceAbundance %>%
  ggplot() +
  geom_violin(aes(x=type, y=medianAbundance), fill="gray80") +
  geom_boxplot(aes(x=type, y=medianAbundance), width=0.1) +
  scale_y_continuous(trans='log10', limits=c(1e-5,1),
                     breaks=trans_breaks('log10', function(x) 10^x),
                     labels=trans_format('log10', math_format(10^.x))) +
  ylab("Species median\nabundance, HMP") +
  DEFAULTS.THEME_PRINT + theme(axis.title.x=element_blank(), axis.text.x=element_text(size=5))
savePNGPDF(paste0(OUTDIR, "8-HMPprevalenceAbundance-residentsColonizers"),
           p8_medianAbundanceResidentColonizers, suppfig8height-0.1, 2.7)
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

# Import species trajectories.
dataSpeciesTrajectories <- read.table("workflow/analysis/annotateSpeciesTrajectoriesSharedStrains/out/speciesAbundances-trajectoriesSharingAnnotated-all.txt.gz",
                                      header=TRUE, stringsAsFactors = FALSE)

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
  p8_species <- dataSpecies %>%
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
    ggtitle(paste0(abxSubject, ", ", gsub("_"," ",substr(ispecies_id, 1, nchar(ispecies_id)-6)))) +
    theme(plot.title=element_text(size=7, margin=margin(0,0,1,0)),
          strip.text.x=element_blank()) +
    theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=4),
          axis.text.y=element_text(size=5)) +
    guides(shape = guide_legend(override.aes = list(size = 1), order=2, reverse=TRUE),
           color = guide_legend(override.aes = list(size = 0.5), order=1)) +
    theme(legend.key.size = unit(0.5, "lines"))
  p8_species
  savePNGPDF(paste0(OUTDIR, "8-", ihh, "-", ispecies_id),
             p8_species + guides(shape="none", color="none"), height, width)
  savePNGPDF(paste0(OUTDIR, "8-speciesTrajectory-legend"),
             get_legend(p8_species), 1.25, 0.5)
}

# Plot examples of strain dynamics
plotSpeciesTrajectory("XI", "Prevotella_copri_61740", 1.4,2.5)
plotSpeciesTrajectory("XK", "Bacteroides_stercoris_56735", 1.4,2.5)


# Supplemental figure 9 - Colonization dynamics on longer timescal --------

suppfig9height <- 1.25

source("workflow/analysis/plotHelpers.R")
OUTDIR <- "workflow/analysis/supp-figures-v2/out/"

# Import the annotated species trajectories.
dataSpeciesTrajectories <- read.table("workflow/analysis/annotateSpeciesTrajectoriesSharedStrains/out/speciesAbundances-trajectoriesSharingAnnotated-all.txt.gz",
                                      header=TRUE, stringsAsFactors = FALSE)
# Plot the cumulative number of colonizing strains during the full study.
# Do not include strains that colonized during the main part of the study.
dataColonizingSpeciesFollowup <- dataSpeciesTrajectories %>%
  mutate(speciesTrajectory=ifelse((!is.na(timeOfColonizationAllTimepoints) | !is.na(timeOfStrainTurnover)),
                                  "colonizedStrainTurnover", "none"),
         speciesTrajectory=ifelse(is.na(speciesTrajectory), "none", speciesTrajectory),
         timeOfStrainColonization=ifelse(!is.na(timeOfColonizationAllTimepoints),
                                         timeOfColonizationAllTimepoints, timeOfStrainTurnover),
         timeOfStrainColonization=ifelse(is.na(timeOfStrainColonization),0,timeOfStrainColonization)) %>%
  dplyr::select(subject, timepoint, species_id, speciesTrajectory, timeOfStrainColonization) %>%
  filter(timepoint>75) %>%
  group_by(subject, timepoint) %>%
  summarize(cumulativeColonization=sum(speciesTrajectory=="colonizedStrainTurnover" & 
                                         timeOfStrainColonization<=timepoint)) %>%
  ungroup() %>% group_by(subject) %>% filter(timepoint==max(timepoint)) %>%
  mutate(subjectResponse=annotateSubjectResponse(subject))
p9_cumulativeNumberColonizingSpecies <- dataColonizingSpeciesFollowup %>%
  mutate(subjectResponseDisplay=gsub(" ","\n",subjectResponse)) %>%
  ggplot() +
  geom_point(aes(x=factor(subjectResponseDisplay, names(PALETTE.ABXRESPONSELINEBREAK)), 
                 y=cumulativeColonization, color=factor(subjectResponse)),
             position=position_jitter(height=0, width=0.2, seed=0), alpha=0.8, size=1) +
  scale_color_manual(values=PALETTE.ABXRESPONSE) + guides(color="none") +
  xlab("Subject response") + ylab("Number of\nnew colonizers") +
  DEFAULTS.THEME_PRINT
savePNGPDF(paste0(OUTDIR, "9-numColonizedTurnover-byResponse-allSamples"), 
           p9_cumulativeNumberColonizingSpecies, suppfig9height, 1.8)
# Test for statistical enrichment of new colonizer events.
# Wilcoxon rank-sum, control vs. minimal/transient: p=0.6107
# Wilcoxon rank-sum, control vs. lasting: p=0.0006968
# Wilcoxon rank-sum, minimal/transient vs. lasting: p=0.001362
wilcox.test(dataColonizingSpeciesFollowup %>% filter(subjectResponse=="control") %>% pull(cumulativeColonization),
            dataColonizingSpeciesFollowup %>% filter(subjectResponse=="transient response") %>% pull(cumulativeColonization))
wilcox.test(dataColonizingSpeciesFollowup %>% filter(subjectResponse=="control") %>% pull(cumulativeColonization),
            dataColonizingSpeciesFollowup %>% filter(subjectResponse=="lasting response") %>% pull(cumulativeColonization))
wilcox.test(dataColonizingSpeciesFollowup %>% filter(subjectResponse=="transient response") %>% pull(cumulativeColonization),
            dataColonizingSpeciesFollowup %>% filter(subjectResponse=="lasting response") %>% pull(cumulativeColonization))

# Import the number of new colonizers detected on longer timescales (i.e. >1 month post-antibiotics).
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
p9_pctColonized_followup <- dataTotalColonizedSummaryEndFollowup %>% 
  mutate(subjectResponseDisplay=gsub(" ","\n",subjectResponse)) %>%
  ggplot() +
  geom_point(aes(x=factor(subjectResponseDisplay, names(PALETTE.ABXRESPONSELINEBREAK)), 
                 y=totalAbundance, color=factor(subjectResponse)),
             position=position_jitter(height=0, width=0.2, seed=0), alpha=0.8, size=1) +
  scale_color_manual(values=PALETTE.ABXRESPONSE) +
  guides(color="none") + ylim(0,1) +
  xlab("Subject response") + ylab("Total abundance,\nnew colonizers") +
  DEFAULTS.THEME_PRINT
savePNGPDF(paste0(OUTDIR,"9-pctColonized-followup"), p9_pctColonized_followup,
           suppfig9height, 1.8)
# Test for statistical enrichment in the total proportion of the gut microbiome 
# that consists of colonizing species.
# Wilcoxon, control vs. minimal/transient: p=0.3562
# Wilcoxon, control vs. lasting: p=0.0009504
# Wilcoxon, transient vs. lasting: p=0.003649
wilcox.test(dataTotalColonizedSummaryEndFollowup %>% filter(subjectResponse=="control") %>% pull(totalAbundance),
            dataTotalColonizedSummaryEndFollowup %>% filter(subjectResponse=="transient response") %>% pull(totalAbundance))
wilcox.test(dataTotalColonizedSummaryEndFollowup %>% filter(subjectResponse=="control") %>% pull(totalAbundance),
            dataTotalColonizedSummaryEndFollowup %>% filter(subjectResponse=="lasting response") %>% pull(totalAbundance))
wilcox.test(dataTotalColonizedSummaryEndFollowup %>% filter(subjectResponse=="transient response") %>% pull(totalAbundance),
            dataTotalColonizedSummaryEndFollowup %>% filter(subjectResponse=="lasting response") %>% pull(totalAbundance))

# Compare the proportion of colonization in subjects with lasting responses
# versus their cohabiting controls.
p9_pctColonized_lastingControls <- dataTotalColonizedSummaryEndFollowup %>%
  filter(substr(subject,1,2) %in% c("XA","XB","XD","XE","XK")) %>%
  mutate(subjectResponseDisplay=gsub(" ","\n",subjectResponse)) %>%
  ggplot() +
  geom_point(aes(x=factor(subjectResponseDisplay, names(PALETTE.ABXRESPONSELINEBREAK)), 
                 y=totalAbundance, color=factor(subjectResponse)),
             position=position_jitter(height=0, width=0.2, seed=0), alpha=0.8, size=1) +
  scale_color_manual(values=PALETTE.ABXRESPONSE) +
  guides(color="none") + ylim(0,1) +
  xlab("Subject") + ylab("Total abundance,\nnew colonizers") +
  DEFAULTS.THEME_PRINT
p9_pctColonized_lastingControls
savePNGPDF(paste0(OUTDIR,"9-pctColonized-lastingResponseControls"), p9_pctColonized_lastingControls,
           suppfig9height, suppfig9height)
# Statistically compare the proportion of colonization in subjects with lasting responses
# and their cohabiting partners.
# Wilcoxon rank-sum test, p=0.0135
wilcox.test(dataTotalColonizedSummaryEndFollowup %>% filter(substr(subject,1,2) %in% c("XA","XB","XD","XE","XK"),
                                                            subjectResponse=="control") %>% pull(totalAbundance),
            dataTotalColonizedSummaryEndFollowup %>% filter(substr(subject,1,2) %in% c("XA","XB","XD","XE","XK"),
                                                            subjectResponse=="lasting response") %>% pull(totalAbundance))


# Plot the transit time in real time.
# Import the trajectory fit summary data.
dataSpeciesTrajectoriesFitSummary <-
  read.table("workflow/analysis/fitColonizationTrajectories/out/speciesTrajectoriesFit-summary.txt",
             header=TRUE, stringsAsFactors = FALSE, sep="\t")
# Plot the distribution of the time in days needed to go
# from undetectable relative abundances to carrying capacity.
PALETTE.NUMTIMEPOINTS <- c("gray","#44AA99","#332288")
names(PALETTE.NUMTIMEPOINTS) <- c("1","2",">2")
p9_transitTime_histogram <- dataSpeciesTrajectoriesFitSummary %>%
  filter(tstarTimepoint>70, typeOfRecoveryColonizationTurnover!="same strain",
         typeOfRecoveryColonizationTurnover!="unknown strain") %>%
  mutate(numTimepointsToReachKcondensed=ifelse(numTimepointsToReachK>=3,">2",as.character(numTimepointsToReachK))) %>%
  mutate(typeOfRecoveryColonizationTurnover=
           ifelse(!is.na(speciesColonizedTransientlyAll) & speciesColonizedTransientlyAll, "new species, transient",
                  ifelse(typeOfRecoveryColonizationTurnover=="new species",
                         "new species, persistent", typeOfRecoveryColonizationTurnover))) %>%
  ggplot() +
  geom_histogram(aes(x=timeIntervalToReachK, 
                     fill=fct_rev(fct_relevel(numTimepointsToReachKcondensed, c("1","2",">2")))),
                 binwidth=50, color="gray30") +
  scale_fill_manual(values=PALETTE.NUMTIMEPOINTS, name="Number of\nintervals") +
  xlab("Transit time (days)") + ylab("Number of\ncolonizers") +
  theme(legend.key.size = unit(0.5, "lines")) +
  guides(fill = guide_legend(override.aes = list(size = 0.5), rev=TRUE), color="none") +
  DEFAULTS.THEME_PRINT + theme(legend.box.margin=margin(0,0,0,-10), legend.position=c(0.7,0.7))
p9_transitTime_histogram
savePNGPDF(paste0(OUTDIR, "9-transitTime-histogram"), p9_transitTime_histogram,
           suppfig9height, 1.5)

# Import the data on colonization cohort clustering.
dataColonizationClustering <-
  read.table("workflow/analysis/analyzeClusteringOfColonizationEvents/out/simulations-colonizationCohorts.txt",
             header=TRUE, stringsAsFactors = FALSE) %>%
  group_by(subject) %>%
  mutate(pvalue=sum(simulatedNumbersOfColonizationEvents>=actualNumberOfColonizationEvents)/n())
# Plot the distribution of simulated colonization clusters
# versus the actual value observed in the data as a violin plot.
p9_colonizationClusteringViolin <- dataColonizationClustering %>%
  ggplot() +
  geom_violin(aes(x=subject, y=simulatedNumbersOfColonizationEvents), fill="grey80", adjust=5) +
  geom_point(aes(x=subject, y=actualNumberOfColonizationEvents,
                 color=factor(pvalue<0.05)), size=1.5, pch=4) +
  scale_color_manual(values=c("black","firebrick3"),
                     name="\n\n\n\n\n\np<0.05?") +
  scale_y_continuous(breaks=seq(0,10,2)) +
  xlab("Subject") + ylab("# of colonization events\ndetected simultaneously") +
  theme(legend.key.size = unit(0.5, "lines")) +
  guides(color = guide_legend(override.aes = list(size = 0.5), reverse=TRUE)) +
  DEFAULTS.THEME_PRINT +
  theme(axis.title=element_text(size=6),
        legend.margin = margin(0, 0, 0, 0, "cm"))
savePNGPDF(paste0(OUTDIR, "9-colonizationClustering-violin"),
           p9_colonizationClusteringViolin, suppfig9height, 2.3)

# Import the XAA B. vulgatus strain frequencies as a portion of the B. vulgatus population.
# This data is calculated from private SNP frequencies.
dataXAABvulgatusStrainFreqs <- read.csv("workflow/analysis/analyzeClusteringOfColonizationEvents/Frequencies_all_strains_from_private_snvs.csv",
                                        header=TRUE, stringsAsFactors = FALSE) %>%
  mutate(sample=gsub("HouseholdTransmission-Stool-","",sample))
# Set strain frequency to 0 at the timepoints with insufficient sequencing coverage.
# This will later be replaced with an unknown strain.
dataXAABvulgatusStrainFreqs <- left_join(samplesRaw %>% filter(subject=="XAA", !(sample %in% sampleBlacklist)),
                                         dataXAABvulgatusStrainFreqs, by=c("sample","timepoint"))
dataXAABvulgatusStrainFreqs <- dataXAABvulgatusStrainFreqs %>%
  complete(sample, strain, fill=list(median_freq_private_snvs=0)) %>%
  filter(!is.na(strain)) %>% mutate(subject="XAA", hh="XA", timepoint=as.integer(substr(sample,5,7)))

# Load the species abundance data.
source("workflow/analysis/generateSpeciesAbundances/loadSpeciesAbundancesFiltered.R")
# Extract the relative abundances of B. vulgatus over time.
dataXAABvulgatusSpeciesAbundances <- speciesAbundances %>%
  filter(subject=="XAA", grepl("Bacteroides_vulgatus", species_id)) %>%
  dplyr::select(subject, timepoint, sample, relative_abundance)
# Impute all of the timepoints at which B. vulgatus was not detected.
dataXAABvulgatusSpeciesAbundances <-
  left_join(samplesRaw %>% filter(subject=="XAA", !(sample %in% sampleBlacklist)),
            dataXAABvulgatusSpeciesAbundances,
            by=c("sample","timepoint","subject")) %>%
  mutate(relative_abundance=ifelse(is.na(relative_abundance), 0, relative_abundance))

# Combine the species and strain abundances of B. vulgatus over time.
dataXAABvulgatusSpeciesAbundances <-
  left_join(dataXAABvulgatusSpeciesAbundances, dataXAABvulgatusStrainFreqs,
            by=c("sample","timepoint","subject","hh"))

# Calculate the relative abundance of each strain by multiplying the strain frequency
# by the relative abundance.
dataXAABvulgatusSpeciesAbundances <- dataXAABvulgatusSpeciesAbundances %>%
  mutate(median_freq_private_snvs=ifelse(is.na(median_freq_private_snvs), 0, median_freq_private_snvs)) %>%
  group_by(timepoint) %>%
  mutate(median_freq_private_snvs_normalized=median_freq_private_snvs/sum(median_freq_private_snvs)) %>%
  mutate(median_freq_private_snvs_normalized=ifelse(is.na(median_freq_private_snvs_normalized), 0,
                                                    median_freq_private_snvs_normalized)) %>%
  mutate(strainAbundance=relative_abundance*median_freq_private_snvs_normalized)

# Where there is insufficient coverage to perform strain-level analysis,
# impute an "unknown" strain that occupies the full relative abundance of the population.
dataXAABvulgatusSpeciesAbundances <- dataXAABvulgatusSpeciesAbundances %>%
  mutate(strainShort=case_when(
    strain=="Major Strain Main Study" ~ "strain1",
    strain=="Minor Strain Main Study" ~ "strain2",
    strain=="Major Strain Follow Up Study" ~ "strain3",
    strain=="Minor Strain Follow Up Study" ~ "strain4"
  )) %>% dplyr::select(subject, timepoint, strainShort, relative_abundance, strainAbundance) %>%
  pivot_wider(names_from=strainShort, values_from=strainAbundance) %>%
  mutate(strain5=relative_abundance-strain1-strain2-strain3-strain4,
         strain5=ifelse(strain5<1e-5,0,strain5)) %>%
  pivot_longer(starts_with("strain"), names_to = "strain", values_to="strainAbundance")

p9_XAABvulgatus <- dataXAABvulgatusSpeciesAbundances %>%
  mutate(timePeriod=ifelse(timepoint<75,"main","follow-up")) %>%
  ggplot() +
  geom_rect(aes(xmin=0, xmax=5, ymin=0, ymax=max(dataXAABvulgatusSpeciesAbundances$relative_abundance)),
            fill=abxColor, alpha=0.8) +
  geom_area(aes(x=timepoint-29, y=strainAbundance, fill=factor(strain))) +
  scale_fill_manual(values=c("#332288","#DDCC77","#882255", "#117733","gray30"), name="") +
  xlab("Study day") + ylab("Relative\nabundance") + guides(fill="none") +
  DEFAULTS.THEME_PRINT
savePNGPDF(paste0(OUTDIR, "9-XAA-Bvulgatus-co-colonization-main"),
           p9_XAABvulgatus + xlim(-30,35), fig4height-0.25, 1.5)
savePNGPDF(paste0(OUTDIR, "9-XAA-Bvulgatus-co-colonization-followup"),
           p9_XAABvulgatus + xlim(35,NA) + DEFAULTS.THEME_NOYAXIS, fig4height-0.25, 1)



# Supplemental figure 10 - XBA dynamics ------------------------------------

source("workflow/analysis/plotHelpers.R")
OUTDIR <- "workflow/analysis/supp-figures-v2/out/"
suppfig10height <- 1.25

# Import species abundances.
source("workflow/analysis/generateSpeciesAbundances/loadSpeciesAbundancesFiltered.R")

# Plot the relative abundance trajectories of the Bacteroides species only
# around the time of the switch between states 2 and 3.
initialBacteroidesXBA <- (speciesAbundances %>% 
                            filter(sample=="XBA-001", grepl("Bacteroides_", species_id),
                                   relative_abundance>1e-4) %>% 
                            arrange(desc(relative_abundance)) %>% dplyr::select(species_id, relative_abundance))$species_id
p10_BacteroidesAbundances <- speciesAbundances %>% group_by(species_id) %>%
  filter(subject=="XBA", max(relative_abundance)>LIMITOFDETECTION,
         !(sample %in% samplesXBAextraFollowup), species_id %in% initialBacteroidesXBA,
         timepoint>300, timepoint<600) %>%
  mutate(relative_abundance=ifelse(relative_abundance<1e-5,1e-5, relative_abundance),
         timepoint=timepoint-29) %>%
  filter(max(relative_abundance)>1e-5) %>%
  ggplot() +
  geom_line(aes(x=factor(timepoint), y=relative_abundance, group=species_id, 
                color=extraShortenSpeciesName(renamePhocaeicola(species_id)))) +
  geom_hline(yintercept=1e-5, linetype="dashed") +
  scale_y_log10() +
  scale_color_manual(values=colorRampPalette(DEFAULTS.PALETTE.TOL.COLORBLINDSAFE)(length(initialBacteroidesXBA)),
                     name="Species") +
  guides(color = guide_legend(override.aes = list(size = 0.3), ncol=3)) +
  theme(legend.key.size = unit(0.3, "lines")) +
  xlab("Study day") + ylab("Relative abundance") +
  DEFAULTS.THEME_PRINT + 
  theme(legend.text=element_text(size=4.5), 
        axis.text.x=element_text(size=4.5, angle=90, hjust=1, vjust=0.5),
        legend.title=element_blank(), legend.position="top")
savePNGPDF(paste0(OUTDIR, "10-abundancesBacteroides"), p10_BacteroidesAbundances + guides(color="none"),
           1.7, 2.4)
savePNGPDF(paste0(OUTDIR, "10-abundancesBacteroides-legend"), get_legend(p10_BacteroidesAbundances),
           0.4, 2)

# Import the full relative abundances of all annotated species.
dataSpeciesTrajectories <-
  read.table("workflow/analysis/annotateSpeciesTrajectoriesSharedStrains/out/speciesAbundances-trajectoriesSharingAnnotated-all.txt.gz",
             header=TRUE, stringsAsFactors = FALSE)
initialBacteroidesXBA <- (dataSpeciesTrajectories %>% 
                            filter(sample=="XBA-001", grepl("Bacteroides_", species_id),
                                   species_id!="Bacteroides_fragilis_54507", relative_abundance>1e-4) %>% 
                            arrange(desc(relative_abundance)) %>% dplyr::select(species_id, relative_abundance))$species_id
XBAspeciesOfInterest <- c("Bacteroides_fragilis_54507", "Eggerthella_lenta_56463",
                          "Bifidobacterium_longum_57796", initialBacteroidesXBA)
XBAspeciesOfInterestShort <- sapply(XBAspeciesOfInterest, shortenSpeciesName)
dataSpeciesTrajectoriesXB <- dataSpeciesTrajectories %>%
  filter(hh=="XB", species_id %in% XBAspeciesOfInterest,
         !(sample %in% samplesXBAextraFollowup),
         (subject=="XBA" | sample %in% c("XBB-001", "XBB-330","XBB-679", "XBB-900"))) %>%
  mutate(species_id_short=gsub("_"," ",substr(species_id, 1, nchar(species_id)-6)))
p10_XBheatmap <- dataSpeciesTrajectoriesXB %>%
  mutate(subject=paste0(subject,"\n"), timepoint=timepoint-29) %>%
  ggplot() +
  geom_tile(aes(x=factor(timepoint), y=fct_relevel(species_id, rev(XBAspeciesOfInterest)), 
                fill=log10(relative_abundance))) +
  geom_point(data=dataSpeciesTrajectoriesXB %>%
               filter(strainTurnover & timepoint>timeOfStrainTurnover),
             aes(x=factor(timepoint), y=fct_relevel(species_id, rev(XBAspeciesOfInterest))),
             color="gray30", shape=1, size=1.5) +
  geom_point(data=dataSpeciesTrajectoriesXB %>%
               mutate(subject=paste0(subject,"\n"), timepoint=timepoint-29) %>%
               group_by(subject, species_id) %>%
               filter(!is.na(sharingHhByTimepoint)) %>%
               filter(sample %in% c(samplesInitial, samplesLastSequencedX)),
             aes(x=factor(timepoint), y=fct_relevel(species_id, rev(XBAspeciesOfInterest)), 
                 shape=factor(sharingHhByTimepoint)), size=0.8) +
  facet_grid(.~subject, scales="free_x", space="free_x") +
  scale_fill_viridis(option="A", name="Log\nrelative\nabundance") +
  scale_color_manual(name="Strain\nturnover?", values=c("black", "#FFFFFF")) +
  scale_shape_manual(values=c(1, 16), name="Strain\nshared?") +
  xlab("Study day") +
  scale_y_discrete(labels=sapply(sapply(XBAspeciesOfInterest, renamePhocaeicola), shortenSpeciesName)) +
  DEFAULTS.THEME_PRINT +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=4),
        axis.text.y=element_text(size=4.5),
        axis.title.y=element_blank()) +
  guides(shape = guide_legend(override.aes = list(size = 0.5), reverse=TRUE)) +
  theme(legend.key.size = unit(0.5, "lines"))
p10_XBheatmap
savePNGPDF(paste0(OUTDIR, "10-XB-heatmap"), p10_XBheatmap, 2, 4.7)



# Import the dynamics of the SNPs that distinguish B. uniformis in XBA and XBB.
dataSNPsBuniformis <- 
  read.table("workflow/analysis/inferXBAstrainDynamicsPersistenceRecolonization/out/Bacteroides_uniformis_57318-diffSitesSubjects-freqs.txt.gz",
             header=TRUE, stringsAsFactors = FALSE)
dataSNPsBuniformisSummary <- dataSNPsBuniformis %>%
  filter(!is.na(polarizedFreq)) %>%
  group_by(subject, timepoint) %>%
  summarize(medianFreq=median(polarizedFreq), pct25Freq=quantile(polarizedFreq, 0.25),
            pct75Freq=quantile(polarizedFreq, 0.75)) 
p10_Buniformis_XBA <- dataSNPsBuniformisSummary %>%
  mutate(species_id="Bacteroides uniformis", timepoint=timepoint-29) %>%
  ggplot() +
  geom_rect(xmin=0, xmax=5, ymin=0, ymax=1, fill=abxColor) +
  geom_ribbon(aes(x=timepoint, ymin=pct25Freq, ymax=pct75Freq, fill=factor(subject)), alpha=0.5) + 
  geom_point(aes(x=timepoint, y=medianFreq, color=factor(subject)), size=0.5, alpha=0.5) +
  geom_line(aes(x=timepoint, y=medianFreq, color=factor(subject))) +
  facet_wrap(~species_id) +
  scale_color_manual(values=c("#882255","#332288"), name="Subject") +
  scale_fill_manual(values=c("#882255","#332288"), name="Subject") +
  xlab("Study day") + ylab("Frequency,\nXBB allele") + ylim(0,1) +
  guides(color = guide_legend(override.aes = list(size = 0.5)),
         fill="none") +
  theme(legend.key.size = unit(0.5, "lines")) +
  DEFAULTS.THEME_PRINT
savePNGPDF(paste0(OUTDIR, "10-Buniformis-XBA"), p10_Buniformis_XBA + xlim(-30,40) + guides(color="none"),
           suppfig10height, 1.5) 
savePNGPDF(paste0(OUTDIR, "10-Buniformis-XBA-followup"), 
           p10_Buniformis_XBA + xlim(90-29,NA) + DEFAULTS.THEME_NOYAXIS,
           suppfig10height, 1.75)


# Import the dynamics of the SNPs that distinguish B. vulgatus in XBA and XBB.
dataSNPsBvulgatus <- 
  read.table("workflow/analysis/inferXBAstrainDynamicsPersistenceRecolonization/out/Bacteroides_vulgatus_57955-diffSitesSubjects-freqs.txt.gz",
             header=TRUE, stringsAsFactors = FALSE)
dataSNPsBvulgatusSummary <- dataSNPsBvulgatus %>%
  filter(!is.na(polarizedFreq)) %>%
  group_by(subject, timepoint) %>%
  summarize(medianFreq=median(polarizedFreq), pct25Freq=quantile(polarizedFreq, 0.25),
            pct75Freq=quantile(polarizedFreq, 0.75)) 
p10_Bvulgatus_XBA <- dataSNPsBvulgatus %>%
  filter(!is.na(polarizedFreq)) %>%
  mutate(species_id="Phocaeicola vulgatus", timepoint=timepoint-29) %>%
  ggplot() +
  geom_rect(xmin=0, xmax=5, ymin=0, ymax=1, fill=abxColor) +
  geom_point(aes(x=timepoint, y=polarizedFreq, color=factor(subject)), size=0.5, alpha=0.5) +
  geom_line(aes(x=timepoint, y=polarizedFreq, color=factor(subject), group=interaction(subject, site_id))) +
  facet_wrap(~species_id) +
  scale_color_manual(values=c("#882255","#332288"), name="Subject") +
  xlab("Study day") + ylab("Frequency,\nXBB allele") + ylim(0,1) +
  guides(color = guide_legend(override.aes = list(size = 0.5))) +
  theme(legend.key.size = unit(0.5, "lines")) +
  DEFAULTS.THEME_PRINT
p10_Bvulgatus_XBA
savePNGPDF(paste0(OUTDIR, "10-Bvulgatus-XBA"), p10_Bvulgatus_XBA + xlim(-30,40) + guides(color="none"),
           suppfig6height, 1.5) 
savePNGPDF(paste0(OUTDIR, "10-Bvulgatus-XBA-followup"), 
           p10_Bvulgatus_XBA + xlim(90-29,NA) + DEFAULTS.THEME_NOYAXIS,
           suppfig6height, 1.75)


# Import the dynamics of the SNPs that distinguish B. vulgatus in XBA and XBB.
dataSNPsBvulgatusSweeps <- 
  read.table("workflow/analysis/inferXBAstrainDynamicsPersistenceRecolonization/out/Bacteroides_vulgatus_57955-diffSitesTimePeriodsXBA-freqs.txt.gz",
             header=TRUE, stringsAsFactors = FALSE)
# Remove sites that have intermediate frequencies in XBB, which may indicate some mapping issue.
p10_Bvulgatus_XBA_sweeps <- dataSNPsBvulgatusSweeps %>%
  filter(!is.na(polarizedFreq)) %>% group_by(site_id) %>%
  filter(min(polarizedFreq[subject=="XBB"])<0.1) %>%
  mutate(species_id="Phocaeicola vulgatus", timepoint=timepoint-29) %>%
  ggplot() +
  geom_rect(xmin=0, xmax=5, ymin=0, ymax=1, fill=abxColor) +
  geom_point(aes(x=timepoint, y=polarizedFreq, color=factor(subject)), size=0.5, alpha=0.5) +
  geom_line(aes(x=timepoint, y=polarizedFreq, color=factor(subject), group=interaction(subject, site_id))) +
  facet_wrap(~species_id) +
  scale_color_manual(values=c("#882255","#332288"), name="Subject") +
  xlab("Study day") + ylab("Frequency,\nXBB allele") + ylim(0,1) +
  guides(color = guide_legend(override.aes = list(size = 0.5))) +
  theme(legend.key.size = unit(0.5, "lines")) +
  DEFAULTS.THEME_PRINT
p10_Bvulgatus_XBA_sweeps
savePNGPDF(paste0(OUTDIR, "10-Bvulgatus-sweeps-XBA-"), p10_Bvulgatus_XBA_sweeps + xlim(-30,40) + guides(color="none"),
           suppfig6height, 1.5) 
savePNGPDF(paste0(OUTDIR, "10-Bvulgatus-sweeps-XBA-followup"), 
           p10_Bvulgatus_XBA_sweeps + xlim(90-29,NA) + DEFAULTS.THEME_NOYAXIS,
           suppfig6height, 1.75)

# Import the dynamics of the SNPs that distinguish B. uniformis in XBA and XBB.
dataSNPsBlongum <- 
  read.table("workflow/analysis/inferXBAstrainDynamicsPersistenceRecolonization/out/Bifidobacterium_longum_57796-diffSitesSubjects-freqs.txt.gz",
             header=TRUE, stringsAsFactors = FALSE)
dataSNPsBlongumSummary <- dataSNPsBlongum %>%
  filter(!is.na(polarizedFreq)) %>%
  group_by(subject, timepoint) %>%
  summarize(medianFreq=median(polarizedFreq), pct25Freq=quantile(polarizedFreq, 0.25),
            pct75Freq=quantile(polarizedFreq, 0.75)) 
p10_Blongum_XBA <- dataSNPsBlongumSummary %>%
  mutate(species_id="B. longum", timepoint=timepoint-29) %>%
  arrange(desc(subject), timepoint) %>%
  ggplot() +
  geom_rect(xmin=0, xmax=5, ymin=0, ymax=1, fill=abxColor) +
  geom_ribbon(aes(x=timepoint, ymin=pct25Freq, ymax=pct75Freq, fill=factor(subject)), alpha=0.5) + 
  geom_point(aes(x=timepoint, y=medianFreq, color=factor(subject)), size=0.5, alpha=0.5) +
  geom_line(aes(x=timepoint, y=medianFreq, color=factor(subject))) +
  facet_wrap(~species_id) +
  scale_color_manual(values=c("#882255","#332288"), name="Subject") +
  scale_fill_manual(values=c("#882255","#332288"), name="Subject") +
  xlab("Study day") + ylab("Frequency,\nXBB allele") + ylim(0,1) +
  guides(color = guide_legend(override.aes = list(size = 0.5)),
         fill="none") +
  theme(legend.key.size = unit(0.5, "lines")) +
  DEFAULTS.THEME_PRINT
p10_Blongum_XBA
savePNGPDF(paste0(OUTDIR, "10-Blongum-XBA"), p10_Blongum_XBA + xlim(-30,40) + guides(color="none"),
           suppfig10height, 1.5) 
savePNGPDF(paste0(OUTDIR, "10-Blongum-XBA-followup"), 
           p10_Blongum_XBA + xlim(90-29,NA) + DEFAULTS.THEME_NOYAXIS,
           suppfig10height, 1.75)

# Plot the correlation between T6SS relative abundance and the relative abundance
# of B. fragilis in XBA.
# Import data on T6SS dynamics.
dataT6SS <- read.table("workflow/analysis/scratch/230305-T6SS/out/T6SScoverage-cleaned.txt",
                       header=TRUE, stringsAsFactors = FALSE)
source("workflow/analysis/generateSpeciesAbundances/loadSpeciesAbundancesFiltered.R")
# Exclude blacklisted samples from further analysis.
dataT6SS <- dataT6SS %>%
  filter((sample %in% samplesRaw$sample), !(sample %in% sampleBlacklist))
# Calculate the number of reads mapping to each bp of the T6SS genes
# per billion reads sequenced.
dataT6SS <- dataT6SS %>%
  mutate(normalizedAbundance=numReads/totalReads/geneLength*1e9)  %>%
  mutate(subject=substr(sample,1,3), timepoint=as.integer(substr(sample,5,7)))
# Extract the data for subjects XBA and XBB.
dataT6SSXB <- dataT6SS %>%
  filter(grepl("XB", sample))
# Extract the data for GA3 effector gene 13 only.
dataT6SSXBA_GA3_E13 <- dataT6SSXB %>%
  filter(GAtype=="GA3", geneNumber==13, subject=="XBA") %>%
  dplyr::select(sample, subject, timepoint, sample, normalizedAbundance)
# Merge the GA3 effector gene normalized abundances
# and the abundances of B. fragilis in XBA at the same timepoints.
dataT6SSXBA_GA3_E13 <- left_join(dataT6SSXBA_GA3_E13,
          speciesAbundances %>% filter(species_id=="Bacteroides_fragilis_54507") %>%
            dplyr::select(sample, subject, timepoint, relative_abundance),
          by=c("sample","subject","timepoint")) %>%
  filter(!is.na(normalizedAbundance), !is.na(relative_abundance))
# Calculate the correlation between B. fragilis relative abundance and
# T6SS G3 normalized abundance.
# Pearson r=0.91, p<2.2e-16
corr.GA3Bfragilis <- cor.test(dataT6SSXBA_GA3_E13 %>% mutate(normalizedAbundance=log10(normalizedAbundance+0.001)) %>% pull(normalizedAbundance),
         dataT6SSXBA_GA3_E13 %>% mutate(relative_abundance=log10(relative_abundance+0.001)) %>% pull(relative_abundance), method="pearson")
# Plot the correlation between B. fragilis relative abundance and GA3 E13 normalized abundance.
p10_T6SSvsBfragilis <- dataT6SSXBA_GA3_E13 %>%
  ggplot() +
  geom_point(aes(x=relative_abundance, y=normalizedAbundance), alpha=0.25, size=1) +
  annotate(geom="text", x=0.002, y=100,
           label=paste0("Pearson r= ", round(corr.GA3Bfragilis$estimate, digits=2),
                        "\np", ifelse(corr.GA3Bfragilis$p.value<1e-16, "< 1e-16",
                                      formatC(corr.GA3Bfragilis$p.value, format="e", digits=2))),
            size=2, hjust=0) +
  scale_x_log10(limits=c(NA, 3)) + scale_y_log10() +
  xlab("Relative abundance,\nBacteroides fragilis") +
  ylab("Normalized abundance,\nGA3-E13") +
  DEFAULTS.THEME_PRINT + theme(axis.title=element_text(size=6))
savePNGPDF(paste0(OUTDIR, "10-Bfragilis-vs-T6SS-GA3-E13"), p10_T6SSvsBfragilis,
           suppfig10height, 1.5)

# Extract the data on T6SS gene abundances in other subjects with lasting abx responses.
dataT6SSlastingResponse <- dataT6SS %>%
  filter(annotateSubjectResponse(subject)=="lasting response",
         subject!="XBA") %>%
  mutate(genotype=paste0(GAtype,"-",geneFunction,geneNumber),
         genotypeShort=paste0(GAtype,"-",geneNumber)) %>%
  group_by(subject, genotype) %>%
  filter(max(normalizedAbundance)>5, geneFunction=="E")
plotT6SS <- function(iSubject, T6SSpalette){
  dataT6SSlastingResponse %>%
    filter(subject==iSubject) %>%
    mutate(timepoint=timepoint-29) %>%
    ggplot() +
    geom_rect(aes(xmin=0, xmax=5, ymin=0, 
              ymax=max(dataT6SSlastingResponse %>% filter(subject==iSubject) %>% 
                         pull(normalizedAbundance))*1.1), fill=abxColor) +
    geom_line(aes(x=timepoint, y=normalizedAbundance, group=genotype,
                  color=factor(genotype))) +
    facet_wrap(~subject) +
    scale_color_manual(values=T6SSpalette, name="Effector") +
    xlab("") + ylab("") + ylim(0,max(dataT6SSlastingResponse %>% filter(subject==iSubject) %>% 
                                       pull(normalizedAbundance))*1.1) +
    guides(color = guide_legend(override.aes = list(size = 0.5)),
           linetype = guide_legend(override.aes = list(size = 0.5))) +
    theme(legend.key.size = unit(0.5, "lines"),
          legend.box.margin = margin(0,0,0,-10)) +
    DEFAULTS.THEME_PRINT 
}
savePNGPDF(paste0(OUTDIR, "10-T6SS-XAA"), 
           plotT6SS("XAA", DEFAULTS.PALETTE.TOL.COLORBLINDSAFE) + 
             guides(color="none") + xlim(-30,35),
           suppfig10height, 1.25) 
savePNGPDF(paste0(OUTDIR, "10-T6SS-XAA-followup"), 
           plotT6SS("XAA", DEFAULTS.PALETTE.TOL.COLORBLINDSAFE) + xlim(35,NA) +
             DEFAULTS.THEME_NOYAXIS,
           suppfig10height, 1.5)
savePNGPDF(paste0(OUTDIR, "10-T6SS-XKA"), 
           plotT6SS("XKA", c(DEFAULTS.PALETTE.TOL.COLORBLINDSAFE[c(4,6,8)])) + 
             guides(color="none") + xlim(-30,35),
           suppfig10height, 1.25) 
savePNGPDF(paste0(OUTDIR, "10-T6SS-XKA-followup"), 
           plotT6SS("XKA", c(DEFAULTS.PALETTE.TOL.COLORBLINDSAFE[c(4,6,8)])) + xlim(35,NA) +
             DEFAULTS.THEME_NOYAXIS,
           suppfig10height, 1.5)
