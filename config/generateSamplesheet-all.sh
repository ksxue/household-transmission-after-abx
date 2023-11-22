# Generate samplesheets for each sample subset.
config/generateSamplesheet-phase1stoolsamples.sh
config/generateSamplesheet-XBstoolsamples-followup.sh
config/generateSamplesheet-phase2stoolsamples.sh
config/generateSamplesheet-phase3stoolsamples.sh
config/generateSamplesheet-phase4stoolsamples.sh

# Merge samplesheets listing raw data files.
cat config/phase1stoolsamples-raw.txt config/XBstoolsamples-followup-raw.txt config/phase2stoolsamples-raw.txt \
  config/phase3stoolsamples-raw.txt config/phase4stoolsamples-raw.txt > config/samples-raw-unfiltered.txt

# Identify the contaminated FASTQ files based on the list of contaminated samples.
config/identifyContaminatedFASTQs.sh

# Extract the FASTQ files that are not listed among the contaminated samples
# from sequencing phases 1, 2, and 3.
grep -F -x -v -f config/contaminatedFASTQs-phase123.txt config/samples-raw-unfiltered.txt \
  > config/samples-raw.txt
