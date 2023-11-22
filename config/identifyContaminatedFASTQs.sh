# Iterate through the list of FASTQ files from phases 1, 2, and 3.
# Iterate through the list of contaminated samples.
# Remove FASTQ files that correspond to the contaminated samples.
rm -f config/contaminatedFASTQs-phase123.txt
while read sample numSuspiciousPairs
do
  grep $sample config/phase1stoolsamples-raw.txt
  grep $sample config/phase2stoolsamples-raw.txt
  grep $sample config/phase3stoolsamples-raw.txt
  grep $sample config/XBstoolsamples-followup-raw.txt
done < <( tail -n +2 config/contaminatedSamples-phase123.txt ) \
  >> config/contaminatedFASTQs-phase123.txt
