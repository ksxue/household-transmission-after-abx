# Generate a list of raw reads for each sample in format sample | read1 | read2.
for f in /oak/stanford/groups/relman/gcloud/compressed_raw_data_backup/dorang/201222-XCBDSalivaMgx/Katherine_Xue/Stool*/*_R1*.fastq.gz
do
  sample=${f##*-Stool-}
  sample=${sample%%_*}
  trimdir="/scratch/groups/relman/kxue/household-transmission-mgx/trimmed/"
  project="HouseholdTransmission-Stool-"
  printf '%s\t%s\t%s\t%s\t%s%s%s-trimmed-pair1.fastq.gz\t%s%s%s-trimmed-pair2.fastq.gz\n' \
    ${sample} ${sample} ${f} ${f/_R1_/_R2_} ${trimdir} ${project} ${sample} ${trimdir} ${project} ${sample}
done > config/XBstoolsamples-followup-raw.txt

