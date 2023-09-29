# This module takes the output of the T6SS read mapping pipeline
# and summarizes the number of reads that map to each gene.
# Summarize reads mapping to each gene.
source("workflow/analysis/summarizeT6SScoverage/generateT6SSCoverageSummary.R")
rm(list=ls())
