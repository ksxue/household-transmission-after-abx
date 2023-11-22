library(tidyverse)
library(data.table)
library(cowplot)

# This script sets up the palettes and functions to plot species relative abundances.

# Import scripts to load species abundance data in the appropriate format for plotting.
source("workflow/analysis/generateSpeciesAbundances/loadSpeciesAbundancesFiltered.R")
# Import the plotting defaults.
source("config/plotDefaults.R")
# Set large color palette.
PALETTE.EXPANDED = colorRampPalette(brewer.pal(11,"Spectral"))
# Set output directory.
OUTDIR <- "workflow/analysis/generateSpeciesAbundances/out/"

# Identify the common species that exceed a relative abundance of 0.5% in any subject
# at some point during the sampling period.
commonSpecies <- speciesAbundances %>%
  group_by(subject, species_id) %>%
  mutate(maxRelativeAbundance=max(relative_abundance)) %>%
  filter(maxRelativeAbundance>0.005)
# Identify the common families that these common species belong to.
commonFamilies <- sort(unique(commonSpecies$family))
commonFamiliesPalette <- c("grey30",PALETTE.EXPANDED(length(commonFamilies)-1))
names(commonFamiliesPalette) <- commonFamilies


# Generate a color palette for common species in this dataset.
# Order the species based on the families that they belong to
# and use the families to assign colors.
commonSpeciesByFamily <- commonSpecies %>% ungroup() %>%
  dplyr::select(species_id, family) %>%
  unique() %>%
  arrange(family, species_id)
commonSpeciesByFamily$color <- commonFamiliesPalette[commonSpeciesByFamily$family]
commonSpeciesByFamily <- commonSpeciesByFamily %>%
  mutate(color=ifelse(family=="","grey30",color))
commonSpeciesList <- commonSpeciesByFamily$species_id
commonSpeciesPalette <- commonSpeciesByFamily$color
names(commonSpeciesPalette) <- commonSpeciesList

# Export the family color palette for reference.
commonFamiliesPalette[["Unknown"]] <- "grey30"
p_palette <- commonSpeciesByFamily %>%
  dplyr::select(family, color) %>% unique() %>%
  ggplot() +
  geom_tile(aes(x=family, y=0, 
                fill=factor(family))) +
  scale_fill_manual(name="Family",
                    labels=c("Unknown",commonFamilies[-1]),
                    values=commonFamiliesPalette) +
  DEFAULTS.THEME_ALL +
  theme(axis.line=element_blank(),
        axis.ticks=element_blank(),
        axis.title=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=10))
save_plot(paste0(OUTDIR,"palette.png"), p_palette, nrow=1.3)

# Export taxonomic lineage of common families.
commonFamilyTaxonomy <- speciesAbundances %>%
  filter(family %in% commonFamilies, family!="") %>%
  dplyr::select(phylum, class, order, family) %>%
  unique() %>%
  dplyr::arrange(phylum, class, order, family)
write.table(commonFamilyTaxonomy, paste0(OUTDIR,"commonFamiliesTaxonomy.txt"),
            row.names=FALSE, quote=FALSE)  

# Export common species and common families palettes.
saveRDS(commonFamiliesPalette, paste0(OUTDIR,"commonFamiliesPalette.rds"))
saveRDS(commonSpeciesPalette, paste0(OUTDIR,"commonSpeciesPalette.rds"))
