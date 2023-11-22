# Use bowtie2 to map reads to a database of T6SS effector, immunity, and structural genes.
# Use the FASTQ files for which reads have already been filtered out.
rule mapReadsT6SS:
	input:
		r1=join(config["filterdir"],"{sample}-filtered.1.fastq.gz"),
		r2=join(config["filterdir"],"{sample}-filtered.2.fastq.gz")
	output:
		bt2log=join(config["T6SSdir"],"{sample}.bt2.log"),
		samfile=join(config["T6SSdir"],"{sample}.sam")
	threads: config['maxCPUs']
	params:
		T6SSref=config['T6SSRef'],
		T6SSdir=config['T6SSdir'],
		project=config['project']
	conda:
		"../../workflow/envs/bioinfobasics-from-history.yml"
	shell:
		"""
		bowtie2 --very-fast-local \
		  -x {params.T6SSref} \
		  -1 {input.r1} -2 {input.r2} \
		  -S {output.samfile} \
		  --al-conc-gz \
		  {params.T6SSdir}/{wildcards.sample}-T6SS-aligned.%.fastq.gz \
		  2> {output.bt2log}
		"""

# Use samtools to sort the SAM files, convert to BAM files (to save space),
# and index the resulting BAM files for use with other tools.
rule sortAndIndexSAMT6SS:
	input:
		samfile=join(config["T6SSdir"],"{sample}.sam")
	output:
		bamfile=join(config["T6SSdir"],"{sample}.bam"),
		idxstats=join(config["T6SSdir"],"{sample}.idxstats")
	threads: config['maxCPUs']
	conda: "../../workflow/envs/bioinfobasics-from-history.yml"
	shell:
		"""
		samtools sort -@ {threads} {input.samfile} > {output.bamfile}
		samtools index {output.bamfile}
		samtools idxstats {output.bamfile} > {output.idxstats}
		rm {input.samfile}
		"""
		
# Use bedtools to calculate sequencing coverage on all T6SS genes for all samples.
rule concatenateCoverageT6SS:
	input:
		expand("workflow/out/T6SS/HouseholdTransmission-Stool-{sample}.bam",sample=samples)
	output:
		covSummary=join(config["T6SSdir"],"T6SS-coverage.txt"),
		sampleList=join(config["T6SSdir"],"T6SS-samples.txt")
	params: T6SSdir=config['T6SSdir']
	conda: "../../workflow/envs/bioinfobasics-from-history.yml"
	shell:
		"""
		ls {params.T6SSdir}/*.bam > {output.sampleList}
		bedtools multicov -q 20 -bams `cat {output.sampleList}` \
		  -bed databases/T6SS-AJV/T6SS_CTD.bed > {output.covSummary}
		"""
		
# Use a shell script to parse the .bt2.log files and extract the number of reads per sample.
# In this case, this refers to the total number of reads in the filtered FASTQ files
# that were used as input for mapping to the T6SS genes.
rule calculateNumReadsT6SS:
	input:
		covSummary=join(config["T6SSdir"],"T6SS-coverage.txt")
	output:
		readSummary=join(config["T6SSdir"],"T6SS-totalReads.txt")
	params: T6SSdir=config['T6SSdir']
	shell:
		"""
		workflow/scripts/countReads.sh {params.T6SSdir} {output.readSummary}
		"""
		