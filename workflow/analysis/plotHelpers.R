library(tidyverse)
library(patchwork)
library(ggrepel)
library(viridis)
library(scales)

# Import plot defaults.
source("workflow/analysis/plotDefaults.R")
# Import study metadata.
source("workflow/analysis/background/background.R")
# Import the sample blacklist.
sampleBlacklist <- read.table("workflow/analysis/background/out/sampleBlacklist.txt",
                              header=TRUE, stringsAsFactors = FALSE)$sample
LIMITOFDETECTION <- 1e-3
# Abx response palette
ABXRESPONSE <- c("control", "minimal response", "transient response", "lasting response")
ABXRESPONSESHORT <- gsub(" response", "", ABXRESPONSE)
#PALETTE.ABXRESPONSE <- c("grey30", "#88CCEE","#44AA99","#CC6677") # based on Tol palette
# Generate a palette with only transient and lasting responses.
PALETTE.ABXRESPONSE <- c("grey30", "#88CCEE","#88CCEE","#CC6677") # based on Tol palette
names(PALETTE.ABXRESPONSE) <- ABXRESPONSE
PALETTE.ABXRESPONSELINEBREAK <- PALETTE.ABXRESPONSE # based on Tol palette
names(PALETTE.ABXRESPONSELINEBREAK) <- gsub(" ", "\n", ABXRESPONSE)
PALETTE.ABXRESPONSESHORT <- PALETTE.ABXRESPONSE
names(PALETTE.ABXRESPONSESHORT) <- ABXRESPONSESHORT
abxColor <- "goldenrod"
ABXCTRLRESPONSE <- c("control", "abx", "altStates")
PALETTE.ABXCONTROLRESPONSE <- c("gray30","goldenrod","coral3")
names(PALETTE.ABXCONTROLRESPONSE) <- ABXCTRLRESPONSE
# Set the colors for species trajectory annotations.
PALETTE.SPECIESTRAJECTORIES <- c("gray30","gray","#117733","#332288","#882255")
names(PALETTE.SPECIESTRAJECTORIES) <- c("notDisrupted","notClassified","recovered","notRecovered","colonized")


# Write a wrapper function that, when given plot parameters,
# exports both a PDF and PNG of the same dimensions using save_plot.
# Height and width are given in inches.
savePNGPDF <- function(name, plot, height, width){
  save_plot(paste0(name,".png"), plot, base_height=height, base_width=width)
  save_plot(paste0(name,".pdf"), plot, base_height=height, base_width=width)
}

# Set theme defaults for plots with axis breaks that do not need
# addition y-axis information.
DEFAULTS.THEME_NOYAXIS <- theme(axis.text.y=element_blank(),
                                axis.title.y=element_blank(),
                                axis.ticks.y=element_blank(),
                                axis.line.y=element_blank())

# Write a function to shorten the species name.
shortenSpeciesName <- function(x){
  gsub("_"," ",substr(x, 1, nchar(x)-6))
}
# Write a function to rename Bacteroides to Phocaeicola where necessary.
renamePhocaeicola <- function(x){
  renamed <- c("Bacteroides_vulgatus_57955", "Bacteroides_massiliensis_44749", "Bacteroides_sartorii_54642")
  ifelse(x %in% renamed, gsub("Bacteroides","Phocaeicola",x), x)
}

# Write a function to shorten MIDAS names to G. species form.
extraShortenSpeciesName <- function(x){
  shortName <- shortenSpeciesName(x)
  paste0(substr(shortName,1,1), ". ", sub(".*? ", "", shortName))
}


# Import the classification of subject responses.
subjectResponses <- read.table("workflow/analysis/classifySubjectResponses/out/subjectResponses.txt",
                               header=TRUE, stringsAsFactors = FALSE, sep="\t")
subjectResponses <- subjectResponses %>%
  mutate(subjectResponse=ifelse(subjectResponse=="minimal response", "transient response", subjectResponse))

# Annotate subject treatment.
annotateSubjectTreatment <- function(x){
  return(ifelse(x %in% subjectsAbx, "abx", "control"))
}
annotateSubjectResponseVector <- subjectResponses$subjectResponse
names(annotateSubjectResponseVector) <- subjectResponses$subject
# Annotate subject response.
annotateSubjectResponse <- function(x){
  return(annotateSubjectResponseVector[x])
}
# Calculate the number of subjects in each response category.
subjectResponsesSummary <- subjectResponses %>%
  group_by(subjectResponse) %>% summarize(numSubjects=n())
subjectResponsesNumSubjects <- subjectResponsesSummary$numSubjects
names(subjectResponsesNumSubjects) <- subjectResponsesSummary$subjectResponse
# Annotate subject response and label number of subjects.
annotateSubjectResponseNumSubjects <- function(x){
  response <- annotateSubjectResponse(x)
  return(paste0(response, " (n=", subjectResponsesNumSubjects[response], ")"))
}
subjectResponseNumSubjectsOrdered <-
  c(annotateSubjectResponseNumSubjects("XBB"), annotateSubjectResponseNumSubjects("XCA"),
    annotateSubjectResponseNumSubjects("XBA"))

# Set a palette for strain and species annotations.
sharingAnnotationPalette <- c("gray30","lightgray","#DDCC77","#882255",
                              "#332288", "#332288") # Based on the Tol color palette.
names(sharingAnnotationPalette) <- c("speciesNotShared","strainsUnknown","strainsNotShared","strainsShared",
                                     "speciesShared-diffHousehold", "speciesShared-sameHousehold")
displayAnnotations <- c("species\nnot shared", "strains\nunknown", "strains\nnot shared", "strains\nshared",
                        "species\nshared,\ndifferent\nhousehold", "species\nshared,\nsame\nhousehold")
names(displayAnnotations) <- names(sharingAnnotationPalette)

# Set a palette for the strain recovery, colonization, and turnover events.
PALETTESTRAINS <- c("gray80","gray30","#44AA99","#AA4499","#AA4499","#AA4499")
names(PALETTESTRAINS) <- c("unknown strain","same strain","new strain","new species",
                           "new species, persistent", "new species, transient")
PALETTESTRAINSFILL <- c("gray80","gray30","#44AA99","#AA4499","#AA4499","white")
names(PALETTESTRAINSFILL) <- names(PALETTESTRAINS)

# Set scale breaks for plots that use days since abx.
timepointScaleBreaks <- c(-15,0,15,30)

# Plot species abundances using days since antibiotics.
# Requires generateSpeciesAbundances/loadSpeciesAbundancePlots.R to be loaded.
plotSpeciesAbundanceDaysFromAbxStart <- function(x){
  x %>%
    ungroup() %>% group_by(subject) %>%
    mutate(timepoint=timepoint-29) %>%
    mutate(minTimepoint=min(timepoint), maxTimepoint=max(timepoint)) %>%
    ggplot() +
    geom_rect(aes(xmin=minTimepoint, xmax=maxTimepoint, ymin=0, ymax=1),
              fill="lightgray") +
    geom_area(aes(x=timepoint, y=relative_abundance,
                  group=fct_relevel(species_id, commonSpecies),
                  fill=factor(color), alpha=alpha),
              color="black", linewidth=0.1) +
    geom_rect(data=x %>% filter(subject %in% subjectsAbx),
              aes(xmin=0, xmax=5), ymin=1,ymax=1.05,
              fill=abxColor) +
    facet_wrap(~subject, ncol=8, scales="free") +
    scale_fill_manual(values=commonSpeciesPalette) +
    scale_alpha_identity() +
    ylim(0,1) + scale_x_continuous(breaks=timepointScaleBreaks) +
    guides(fill="none") +
    xlab("Study day") + ylab("Relative abundance") +
    DEFAULTS.THEME_ALL
}
