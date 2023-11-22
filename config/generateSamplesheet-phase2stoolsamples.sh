# Generate a samplesheet in which each sample is listed on a separate line,
# with separate lines for different sequencing lanes in format sample | read1 | read2.
# First do this for the first lane of sequencing data for the phase 2 stool samples.
for f in /oak/stanford/groups/relman/gcloud/compressed_raw_data_backup/kxue/210402-Phase2MgxStool-Lane1/*Stool*_R1_*.fastq.gz
do
  sample=${f##*/}
  sample=${sample:28:7}
  sample=${sample/_/-}
  sample=${sample/_/-}
  sample=${sample/_/-}
  samplelane="${sample}-L001"
  trimdir="/scratch/groups/relman/kxue/household-transmission-mgx/trimmed/"
  project="HouseholdTransmission-Stool-"
  printf '%s\t%s\t%s\t%s\t%s%s%s-trimmed-pair1.fastq.gz\t%s%s%s-trimmed-pair2.fastq.gz\n' \
    ${sample} ${samplelane} ${f} ${f/_R1_/_R2_} ${trimdir} ${project} ${samplelane} ${trimdir} ${project} ${samplelane}
done > config/phase2stoolsamples-raw.tmp

# Also do this for the second lane of sequencing data for the phase 2 stool samples.
for f in /oak/stanford/groups/relman/gcloud/compressed_raw_data_backup/kxue/210415-Phase2MgxStool-Lane2/*Stool*_R1_*.fastq.gz
do
  sample=${f##*/}
  sample=${sample:28:7}
  sample=${sample/_/-}
  sample=${sample/_/-}
  sample=${sample/_/-}
  samplelane="${sample}-L002"
  trimdir="/scratch/groups/relman/kxue/household-transmission-mgx/trimmed/"
  project="HouseholdTransmission-Stool-"
  printf '%s\t%s\t%s\t%s\t%s%s%s-trimmed-pair1.fastq.gz\t%s%s%s-trimmed-pair2.fastq.gz\n' \
    ${sample} ${samplelane} ${f} ${f/_R1_/_R2_} ${trimdir} ${project} ${samplelane} ${trimdir} ${project} ${samplelane}
done >> config/phase2stoolsamples-raw.tmp

# Also do this for the third lane of sequencing data for the phase 2 stool samples.
for f in /oak/stanford/groups/relman/gcloud/compressed_raw_data_backup/kxue/210602-Phase2MgxStool-Lane3/*Stool*_R1_*.fastq.gz
do
  sample=${f##*/}
  sample=${sample:28:7}
  sample=${sample/_/-}
  sample=${sample/_/-}
  sample=${sample/_/-}
  samplelane="${sample}-L003"
  trimdir="/scratch/groups/relman/kxue/household-transmission-mgx/trimmed/"
  project="HouseholdTransmission-Stool-"
  printf '%s\t%s\t%s\t%s\t%s%s%s-trimmed-pair1.fastq.gz\t%s%s%s-trimmed-pair2.fastq.gz\n' \
    ${sample} ${samplelane} ${f} ${f/_R1_/_R2_} ${trimdir} ${project} ${samplelane} ${trimdir} ${project} ${samplelane}
done >> config/phase2stoolsamples-raw.tmp

# Remove the test samples and XNB-064, which produced zero reads.
grep -v test config/phase2stoolsamples-raw.tmp | grep -v XNB_064 | sort > config/phase2stoolsamples-raw.txt
rm -f config/phase2stoolsamples-raw.tmp
