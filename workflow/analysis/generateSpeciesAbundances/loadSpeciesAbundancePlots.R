# This script sets up the plotting environment and plotting functions
# to generate species abundance plots.

library(tidyverse)
library(data.table)


# Import plot defaults.
source("config/plotDefaults.R")
# Import background info on household cohort.
source("workflow/analysis/background/background.R")
# Import species abundances.
source("workflow/analysis/generateSpeciesAbundances/loadSpeciesAbundancesFiltered.R")
# Import common species palette.
commonSpeciesPalette <- readRDS("workflow/analysis/generateSpeciesAbundances/out/commonSpeciesPalette.rds")
commonSpecies <- names(commonSpeciesPalette)
# Import common family palette.
commonFamiliesPalette <- readRDS("workflow/analysis/generateSpeciesAbundances/out/commonFamiliesPalette.rds")

# Add standard color and alpha variables for plotting.
speciesAbundances <- speciesAbundances %>%
  mutate(color=species_id, alpha=1)

# Pare down the species abundance matrix to include only common species.
speciesAbundancesCommon <- speciesAbundances %>%
  filter(species_id %in% commonSpecies)

plotSpeciesAbundance <- function(x){
  x %>%
    ungroup() %>% group_by(subject) %>%
    mutate(minTimepoint=min(timepoint), maxTimepoint=max(timepoint)) %>%
    ggplot() +
    geom_rect(aes(xmin=minTimepoint, xmax=maxTimepoint, ymin=0, ymax=1),
              fill="lightgray") +
    geom_area(aes(x=timepoint, y=relative_abundance,
                  group=fct_relevel(species_id, commonSpecies),
                  fill=factor(color), alpha=alpha),
              color="black", linewidth=0.1) +
    geom_rect(data=x %>% filter(subject %in% subjectsAbx),
              aes(xmin=29, xmax=34), ymin=1,ymax=1.05, fill="goldenrod") +
    facet_wrap(~subject, ncol=8, scales="free") +
    scale_fill_manual(values=commonSpeciesPalette) +
    scale_alpha_identity() +
    ylim(0,1) +
    guides(fill="none") +
    xlab("Study day") + ylab("Relative abundance") +
    DEFAULTS.THEME_ALL
}
