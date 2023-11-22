library(tidyverse)

# Set output directory.
OUTDIR <- "workflow/analysis/calculateHMPspeciesPrevalenceAbundance/out/"

# Import the list of HMP species prevalences and abundances.
speciesPrevalenceAbundanceHMP <- 
  read.table("workflow/analysis/calculateHMPspeciesPrevalenceAbundance/out/HMP-speciesPrevalenceAbundance.txt",
             header=TRUE, stringsAsFactors = FALSE)
# Import the list of resident species.
residentSpecies <-
  read.table("workflow/analysis/annotateSpeciesTrajectoriesSharedStrains/out/residentSpecies.txt",
             header=TRUE, stringsAsFactors = FALSE)
# Import the list of colonizing strains and species.
colonizingStrainsSpecies <-
  read.table("workflow/analysis/annotateSpeciesTrajectoriesSharedStrains/out/colonizingStrainsSpecies.txt",
             header=TRUE, stringsAsFactors = FALSE)

# Attach the HMP prevalence and abundance to the resident and colonizing species.
residentSpeciesAnnotated <- 
  left_join(residentSpecies, speciesPrevalenceAbundanceHMP,
            by=c("species_id"))
colonizingSpeciesAnnotated <- 
  left_join(colonizingStrainsSpecies, speciesPrevalenceAbundanceHMP,
            by=c("species_id"))

# Combine the distributions of resident and colonizing species.
residentColonizingSpeciesAnnotated <-
  rbind(residentSpeciesAnnotated %>% dplyr::select(species_id, prevalence, medianAbundance) %>%
          mutate(type="resident\nspecies"),
        colonizingSpeciesAnnotated %>% 
          mutate(type=ifelse(timeOfColonization<75, "new colonizer\n<1 month post-abx", "new colonizer\n>1 month post-abx")) %>%
          dplyr::select(species_id, prevalence, medianAbundance, type))
write.table(residentColonizingSpeciesAnnotated, 
            paste0(OUTDIR, "residentColonizingSpecies-prevalenceAbundance.txt"),
            row.names=FALSE, quote=TRUE)
