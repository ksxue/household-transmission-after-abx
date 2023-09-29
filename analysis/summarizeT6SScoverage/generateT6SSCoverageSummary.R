library(tidyverse)

# Set working directory.
dir <- "workflow/analysis/summarizeT6SScoverage/out/"

# Import list of samples analyzed.
dataSamples <- read.table("workflow/out/T6SS/T6SS-samples.txt",
                      header=FALSE, stringsAsFactors = FALSE)
colnames(dataSamples) <- c("sample")
# Trim off the file path and file extension information.
samples <- (dataSamples %>%
  mutate(sampleshort=gsub(".*/","",sample),
         sampleshort=gsub(".bam","",sampleshort)))$sampleshort

# Import total read coverage on mapped samples.
dataTotalReads <- read.table("workflow/out/T6SS/T6SS-totalReads.txt",
                             header=FALSE, stringsAsFactors = FALSE)
colnames(dataTotalReads) <- c("sample","totalReads")

# Import coverage summary for T6SS genes.
dataT6SScoverage <- read.table("workflow/out/T6SS/T6SS-coverage.txt",
                               header=FALSE, stringsAsFactors = FALSE)
colnames(dataT6SScoverage) <- c("chrom","chromStart","chromEnd",samples)


# Tidy the coverage summary and add info on total read coverage.
dataT6SS <- dataT6SScoverage %>%
  pivot_longer(all_of(samples), names_to="sample", values_to="numReads") %>%
  left_join(dataTotalReads, by=c("sample"))
# Clean up the dataframe by removing unnecessary columns and 
# renaming remaining columns to be more clear.
# Also simplify the sample names.
dataT6SS <- dataT6SS %>%
  dplyr::rename(gene=chrom, geneLength=chromEnd) %>%
  mutate(sample=gsub("HouseholdTransmission-Stool-","",sample)) %>%
  dplyr::select(sample, gene, geneLength, numReads, totalReads)
# Parse the gene names to identify the GA system,
# classify genes as effector, immunity, or structural genes,
# and identify the gene name.
dataT6SS <- dataT6SS %>%
  mutate(GAtype=substr(gene,1,3),
         geneFunction=ifelse(grepl("_E",gene),"E",
                  ifelse(grepl("_I",gene),"I","S")),
         geneNumber=ifelse(geneFunction=="E",gsub("_.*","", gsub(".*_E","",gene)),
                    ifelse(geneFunction=="I",gsub(".*_I","",gene),
                           gsub(".*_","",gene))))
# Remove underscores from the gene number.
dataT6SS <- dataT6SS %>%
  mutate(geneNumber=gsub("_", "", geneNumber))

# Check that the GA type, gene type, and gene number are identified correctly.
dataT6SSGenes <- dataT6SS %>%
  dplyr::select(gene, GAtype, geneFunction, geneNumber) %>%
  unique()

# Export the dataframe corresponding to the T6SS coverage.
write.table(dataT6SS %>% 
              dplyr::select(sample, gene, GAtype, geneFunction, geneNumber,
                            geneLength, numReads, totalReads), 
            paste0(dir,"T6SScoverage-cleaned.txt"),
            row.names=FALSE, quote=FALSE)
