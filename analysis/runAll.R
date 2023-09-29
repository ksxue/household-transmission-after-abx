# Run all of the modules in the data analysis pipeline.

# MODULE 0. SET UP ANALYSIS AND GENERATE LIST OF BLACKLISTED SAMPLES.

# First, identify contaminated samples.
# Doing this requires the output of the calculateFixedDiffs snakemake rule.
source("workflow/analysis/identifyContaminatedSamples/runAll-identifyContaminatedSamples.R")
# Next, import the background information for the sequenced samples
# and generate the list of blacklisted samples.
source("workflow/analysis/background/runAll-background.R")
# Parse participant metadata for downstream analyses.
source("workflow/analysis/parseParticipantMetadata/runAll-parseParticipantMetadata.R")


# MODULE 1. GENERATE SPECIES ABUNDANCES AND ANALYZE COMMUNITY DIVERSITY.

# Generate the dataframe of species abundances.
source("workflow/analysis/generateSpeciesAbundances/runAll-generateSpeciesAbundances.R")
# Compare the species abundances of samples collected at the initial sampling timepoint.
source("workflow/analysis/compareSpeciesAbundances/runAll-compareSpeciesAbundances.R")
# Calculate species alpha, beta, and gamma diversity.
source("workflow/analysis/calculateDiversityStats/runAll-calculateDiversityStats.R")
# Calculate subject responses based on beta diversity.
source("workflow/analysis/classifySubjectResponses/runAll-classifySubjectResponses.R")


# MODULE 2. INFER SHARED STRAINS AND CONSOLIDATE CALLS.
# Infer shared strains using the fixed differences method.
source("workflow/analysis/identifySharedStrains-fixedDiffs/runAll-identifySharedStrains-fixedDiffs.R")
# Infer shared strains using the strain fishing method.
source("workflow/analysis/identifySharedStrains-strainFishing/runAll-identifySharedStrains-strainFishing.R")
# Consolidate the calls of shared strains from the two methods.
source("workflow/analysis/compareSharedStrainCalls/runAll-compareSharedStrainCalls.R")

# MODULE 3. ANNOTATE SPECIES TRAJECTORIES AND SHARED STRAINS
# Directly compare the abundances of shared species and strains
# at certain timepoints of interest.
source("workflow/analysis/integrateSpeciesAbundancesSharedStrains/runAll-integrateSpeciesAbundancesSharedStrains.R")
# Annotate species trajectories and strain sharing, both within a subject
# and between subjects.
source("workflow/analysis/annotateSpeciesTrajectoriesSharedStrains/runAll-annotateSpeciesTrajectoriesSharedStrains.R")

# MODULE 4. ANALYZE PROPERTIES OF SPECIES TRAJECTORIES AND STRAIN SHARING.
# Calculate the taxonomic enrichment of disrupted, recovered, and colonizing species,
# along with shared strains.
source("workflow/analysis/calculateTaxonomicEnrichment/runAll-calculateTaxonomicEnrichment.R")
# Calculate the prevalence and abundance of species in the HMP.
# Also calculate the prevalence and abundance of resident and colonizing species in the HMP.
source("workflow/analysis/calculateHMPspeciesPrevalenceAbundance/runAll-calculateHMPspeciesPrevalenceAbundance.R")

# MODULE 5. FIT TRAJECTORIES TO A SIMPLE MODEL OF EXPONENTIAL GROWTH
# TO INFER CARRYING CAPACITY (K), SATURATION TIME (t*), COLONIZATION DELAY,
# AND TIMING OF COLONIZATION.
# Parse the trajectory fits to infer key parameters.
source("workflow/analysis/fitColonizationTrajectories/runAll-fitColonizationTrajectories.R")
# Analyze the timing of colonization events using permutation tests.
source("workflow/analysis/analyzeClusteringOfColonizationEvents/runAll-analyzeClusteringOfColonizationEvents.R")


# MODULE 6. INFER THE STRAIN DYNAMICS IN XBA.
# Infer strain dynamics in XBA to distinguish between persistence and recolonization.
source("workflow/analysis/inferXBAstrainDynamicsPersistenceRecolonization/runAll-inferXBAstrainDynamicsPersistenceRecolonization.R")

# MODULE 7. ANALYZE THE DYNAMICS OF T6SS.
# Summarize the reads mapping to T6SS genes.
source("workflow/analysis/summarizeT6SScoverage/runAll-generateT6SSCoverageSummary.R")

# MODULE 8. COMPARE THE ABUNDANCE AND PREVALENCE OF RESIDENT SPECIES AND NEW COLONIZERS
# IN THE HUMAN MICROBIOME PROJECT.
# Summarize the prevalence and abundance of species in the HMP and compare to 
# resident species and new colonizers.
source("workflow/analysis/calculateHMPspeciesPrevalenceAbundance/runAll-calculateHMPspeciesPrevalenceAbundance.R")

