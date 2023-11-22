# Generate a samplesheet in which each sample is listed on a separate line,
# with separate lines for different sequencing lanes in format sample | read1 | read2.
# First do this for the first lane of sequencing data for the phase 3 stool samples.
for f in /oak/stanford/groups/relman/gcloud/compressed_raw_data_backup/kxue/210603-Phase3MgxStool-Lane1/*Stool*_R1_*.fastq.gz
do
  sample=${f##*/}
  sample=${sample:28:7}
  sample=${sample/_/-}
  sample=${sample/_/-}
  sample=${sample/_/-}
  samplelane="${sample}-L006"
  trimdir="/scratch/groups/relman/kxue/household-transmission-mgx/trimmed/"
  project="HouseholdTransmission-Stool-"
  printf '%s\t%s\t%s\t%s\t%s%s%s-trimmed-pair1.fastq.gz\t%s%s%s-trimmed-pair2.fastq.gz\n' \
    ${sample} ${samplelane} ${f} ${f/_R1_/_R2_} ${trimdir} ${project} ${samplelane} ${trimdir} ${project} ${samplelane}
done > config/phase3stoolsamples-raw.tmp

# Also do this for the second lane of sequencing data for the phase 3 stool samples.
# Note that this second lane contains separate samples compared to the first lane,
# though they will still be annotated with the L002 designation for clarity.
for f in /oak/stanford/groups/relman/gcloud/compressed_raw_data_backup/kxue/220127-mgx-Phase3MgxStool-Followup/*Stool*_R1_*.fastq.gz
do
  sample=${f##*/}
  sample=${sample:28:7}
  sample=${sample/_/-}
  sample=${sample/_/-}
  sample=${sample/_/-}
  samplelane="${sample}-L007"
  trimdir="/scratch/groups/relman/kxue/household-transmission-mgx/trimmed/"
  project="HouseholdTransmission-Stool-"
  printf '%s\t%s\t%s\t%s\t%s%s%s-trimmed-pair1.fastq.gz\t%s%s%s-trimmed-pair2.fastq.gz\n' \
    ${sample} ${samplelane} ${f} ${f/_R1_/_R2_} ${trimdir} ${project} ${samplelane} ${trimdir} ${project} ${samplelane}
done >> config/phase3stoolsamples-raw.tmp

# Remove the blank wells from downstream analysis.
# Also remove files that contained no sequencing reads, i.e. had the L00M designation.
grep -v BLANK config/phase3stoolsamples-raw.tmp | grep -v XBA-038 | grep -v XKA-027 | grep -v XDB_585 | \
  sort > config/phase3stoolsamples-raw.txt
rm -f config/phase3stoolsamples-raw.tmp
