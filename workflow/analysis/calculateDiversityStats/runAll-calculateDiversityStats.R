# Import the filtered species abundances, generate a phyloseq object,
# and calculate alpha, beta, and gamma diversity statistics on the phyloseq.

# Load the filtered species abundances and output phyloseq objects.
# Output one phyloseq object with all of the observations listed
# and a second that filters for only observations above a relative abundance of 1e-3.
source("workflow/analysis/calculateDiversityStats/generateSpeciesPhyloseq.R")
rm(list=ls())

# Load the filtered species abundances and output family-level phyloseq objects.
# Output one phyloseq object with all of the observations listed
# and a second that filters for only observations above a relative abundance of 1e-3.
source("workflow/analysis/calculateDiversityStats/generateFamilyPhyloseq.R")
rm(list=ls())

# Calculate the alpha diversity on the phyloseq object with and without
# filtering on species abundances.
source("workflow/analysis/calculateDiversityStats/generateAlphaDiversity.R")
rm(list=ls())

# Calculate the median pre-antibiotic alpha diversity based on various statistics
# and the time required to recover that diversity after antibiotics.
source("workflow/analysis/calculateDiversityStats/calculateTimeToRecoverDiversity.R")
rm(list=ls())

# Calculate the beta diversity on the phyloseq object using various statistics.
source("workflow/analysis/calculateDiversityStats/generateBetaDiversity.R")
rm(list=ls())

# Calculate the compositional variability on the phyloseq object.
source("workflow/analysis/calculateDiversityStats/generateCompositionalVariability.R")
rm(list=ls())
