# Generate a samplesheet in which each sample is listed on a separate line,
# with separate lines for different sequencing lanes in format sample | read1 | read2.
# First do this for the first lane of sequencing data for the phase 3 stool samples.
for f in /oak/stanford/groups/relman/gcloud/compressed_raw_data_backup/kxue/221121-mgx-Phase4MgxStool-Followup/*Stool*_R1_*.fastq.gz
do
  sample=${f##*/}
  sample=${sample:28:7}
  sample=${sample/_/-}
  sample=${sample/_/-}
  sample=${sample/_/-}
  samplelane="${sample}-phase4"
  replicate=${f##*/}
  replicate=${replicate:36:1}
  if [[ ${replicate} == "2" ]]
  then
    samplelane="${sample}-phase4-2"
  fi
  trimdir="/scratch/groups/relman/kxue/household-transmission-mgx/trimmed/"
  project="HouseholdTransmission-Stool-"
  printf '%s\t%s\t%s\t%s\t%s%s%s-trimmed-pair1.fastq.gz\t%s%s%s-trimmed-pair2.fastq.gz\n' \
    ${sample} ${samplelane} ${f} ${f/_R1_/_R2_} ${trimdir} ${project} ${samplelane} ${trimdir} ${project} ${samplelane}
done > config/phase4stoolsamples-raw.tmp

# Remove the blank wells from downstream analysis.
# Also remove files that contained no sequencing reads, i.e. had the L00M designation.
grep -v BLANK config/phase4stoolsamples-raw.tmp | grep -v L00M | \
  sort > config/phase4stoolsamples-raw.txt
rm -f config/phase4stoolsamples-raw.tmp
