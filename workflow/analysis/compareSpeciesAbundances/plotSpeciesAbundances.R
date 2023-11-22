# This script plots the species abundances for the main study and follow-up sampling
# for each subject.

library(tidyverse)

# Import the study metadata.
source("workflow/analysis/background/background.R")
# Import the plotting defaults.
source("workflow/analysis/plotDefaults.R")

# Set the output directory.
OUTDIR <- "workflow/analysis/generateSpeciesAbundances/out/speciesAbundancePlots/"

# Import the species abundances and the plotting defaults.
source("workflow/analysis/generateSpeciesAbundances/loadSpeciesAbundancePlots.R")

# For each subject in study arm X, plot species abundances during the main part of the study.
foreach(isubject=subjectsX, .combine="c") %do% {
  p <- plotSpeciesAbundance(speciesAbundances %>% 
                              filter(subject==isubject, sample %in% samplesXmain))
  save_plot(paste0(OUTDIR,isubject,"-main.png"), p, ncol=0.7, nrow=0.7)
}

# For each subject in study arm X, plot species abundances during the full sampling period,
# including follow-up sampling.
foreach(isubject=subjectsX, .combine="c") %do% {
  p <- plotSpeciesAbundance(speciesAbundances %>% 
                              filter(subject==isubject))
  numSamplingDays <- max((speciesAbundances %>% 
                            filter(subject==isubject))$timepoint)
  # Scale the width of the plot based on the number of sampling days.
  save_plot(paste0(OUTDIR,isubject,"-followup.png"), p, 
            ncol=0.3*numSamplingDays/64, nrow=0.7)
}

# Plot a summary of all of the species abundances during the main study
# for subjects in study arm X.
p <- plotSpeciesAbundance(speciesAbundances %>%
                            filter(sample %in% samplesXmain)) +
  facet_wrap(~subject, ncol=8, scales="free")
save_plot(paste0(OUTDIR,"studyArmX-main.png"), p, 
          ncol=0.7*8, nrow=0.7*6)
# Plot a summary of all of the species abundances during the main study
# for antibiotic-taking subjects.
p <- plotSpeciesAbundance(speciesAbundances %>%
                            filter(sample %in% samplesXmain,
                                   subject %in% subjectsAbx)) +
  facet_wrap(~subject, ncol=6, scales="free")
save_plot(paste0(OUTDIR,"studyArmX-main-abx.png"), p, 
          ncol=0.7*6, nrow=0.7*4)

# For each subject in study arms Y and Z, plot species abundances during the full sampling period,
# including follow-up sampling.
foreach(isubject=c(subjectsY, subjectsZ), .combine="c") %do% {
  p <- plotSpeciesAbundance(speciesAbundances %>% 
                              filter(subject==isubject))
  numSamplingDays <- max((speciesAbundances %>% 
                            filter(subject==isubject))$timepoint)
  # Scale the width of the plot based on the number of sampling days.
  save_plot(paste0(OUTDIR,isubject,"-all.png"), p, 
            ncol=max(0.7,0.3*numSamplingDays/64), nrow=0.7)
}
