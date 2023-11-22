# Generate a samplesheet in which each sample is listed on a separate line,
# with separate lines for different sequencing lanes in format sample | read1 | read2.
for f in /oak/stanford/groups/relman/gcloud/compressed_raw_data_backup/kxue/200917-Phase1MgxStool/*Stool*_R1.fastq.gz
do
  samplelane=${f##*-Stool-}
  samplelane=${samplelane%%_R1*}
  sample=${samplelane%%_*}
  trimdir="/scratch/groups/relman/kxue/household-transmission-mgx/trimmed/"
  project="HouseholdTransmission-Stool-"
  printf '%s\t%s\t%s\t%s\t%s%s%s-trimmed-pair1.fastq.gz\t%s%s%s-trimmed-pair2.fastq.gz\n' \
    ${sample} ${samplelane} ${f} ${f/_R1/_R2} ${trimdir} ${project} ${samplelane} ${trimdir} ${project} ${samplelane}
done > config/phase1stoolsamples-raw.tmp

# Remove the test samples.
grep -v test config/phase1stoolsamples-raw.tmp > config/phase1stoolsamples-raw.txt
rm -f config/phase1stoolsamples-raw.tmp
