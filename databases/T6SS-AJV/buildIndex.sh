# This script is meant to be run from the top level of the Github repo.
# It builds a bowtie2 index from FASTA files containing T6SS genes.
# This script assumes that the conda environment bioinfobasics is loaded.

dir="databases/T6SS-AJV/"

# Concatenate the sequences of the CTD of the effectors, the immunity genes,
# and the structural genes.
cat ${dir}/*_effectors_CTD.fasta ${dir}/*_immunity.fasta ${dir}/*_structural.fasta \
  > ${dir}/T6SS_CTD.fasta
  
# Build a bowtie2 index from the T6SS sequences.
bowtie2-build -f ${dir}/T6SS_CTD.fasta ${dir}/T6SS_CTD
