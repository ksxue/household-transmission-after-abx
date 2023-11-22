# Use skewer to trim adapters from raw reads.
rule trimAdapters:
	input:
		r1= lambda wildcards: df.query('samplelane == @wildcards.sample')['read1'].tolist()[0],
		r2= lambda wildcards: df.query('samplelane == @wildcards.sample')['read2'].tolist()[0],
	params:
		outdir=config['trimdir'],
		project=config['project']
	output:
		trim1=join(config["trimdir"],"".join([config["project"],"-{sample}-trimmed-pair1.fastq.gz"])),
		trim2=join(config["trimdir"],"".join([config["project"],"-{sample}-trimmed-pair2.fastq.gz"])),
		done=join(config["trimdir"],"".join([config["project"],"-{sample}.done"]))
	threads: config['maxCPUs']
	conda:
		"../../workflow/envs/skewer-no-builds.yml"
	shell:
		"""
		skewer -x CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -y CTGTCTCTTATACACATCTGACGCTGCCGACGA \
			-t {threads} -z -o {params.outdir}/{params.project}-{wildcards.sample} {input.r1} {input.r2}
		touch {output.done}
		"""

# Use bowtie2 to filter out reads that map to the human genome.
# At this stage, also combine sequencing reads from all lanes.
rule filterOutHumanReads:
	input:
		trim1= lambda wildcards: dict(tuple(df.groupby(['sample'])))[wildcards.sample]['trim1'].tolist(),
		trim2= lambda wildcards: dict(tuple(df.groupby(['sample'])))[wildcards.sample]['trim2'].tolist()
	output:
		filtered1=join(config["filterdir"],"".join([config['project'],"-{sample}-filtered.1.fastq.gz"])),
		filtered2=join(config["filterdir"],"".join([config['project'],"-{sample}-filtered.2.fastq.gz"])),
		bt2log=join(config["filterdir"],"".join([config['project'],"-{sample}.bt2.log"])),
		samfile=temp(join(config["filterdir"],"".join([config['project'],"-{sample}.sam"])))
	threads: config['maxCPUs']
	params:
		humanref=config['humanGenomeRef'],
		filterdir=config['filterdir'],
		project=config['project'],
		trim1= lambda wildcards: ",".join(dict(tuple(df.groupby(['sample'])))[wildcards.sample]['trim1'].tolist()),
		trim2= lambda wildcards: ",".join(dict(tuple(df.groupby(['sample'])))[wildcards.sample]['trim2'].tolist())
	conda:
		"../../workflow/envs/bowtie2-no-builds.yml"
	shell:
		"""
		bowtie2 --very-fast \
		  -x {params.humanref} \
		  -1 {params.trim1} -2 {params.trim2} \
		  -S {params.filterdir}/{params.project}-{wildcards.sample}.sam \
		  --un-conc-gz {params.filterdir}/{params.project}-{wildcards.sample}-filtered \
		  2> {output.bt2log}
		mv {params.filterdir}/{params.project}-{wildcards.sample}-filtered.1 {params.filterdir}/{params.project}-{wildcards.sample}-filtered.1.fastq.gz
		mv {params.filterdir}/{params.project}-{wildcards.sample}-filtered.2 {params.filterdir}/{params.project}-{wildcards.sample}-filtered.2.fastq.gz
		"""
