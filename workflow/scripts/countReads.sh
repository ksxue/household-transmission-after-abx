# Given a directory with bowtie2 output files with the extension .bt2.log,
# this script parses the .bt2.log files to extract the number of reads in the original file.
# Note that this is not the number of reads mapped, but the total number of reads.

dir=$1
output=$2

for f in $1/*.bt2.log
do
  sample=${f##*/}
  sample=${sample%%.*}
  echo ${sample} $(grep 'reads;' ${f} | cut -f1 -d' ')
done > $2
