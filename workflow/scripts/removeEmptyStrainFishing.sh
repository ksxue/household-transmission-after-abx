# This script iterates through the strain fishing output.
# It removes empty gzipped output files to prevent them from interfering with downstream operations.

for f in workflow/report/performStrainFishing/*/strainFishing.txt.gz
do
  size=$(zless -S $f | wc -l )
  if [ $size -eq 1 ]
  then
	echo "Remove $f"
	rm $f
  fi
done

touch workflow/report/performStrainFishing/removeEmptyStrainFishing.done