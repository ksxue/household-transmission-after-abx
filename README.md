# household-transmission-after-abx

The code in this repository performs the analyses described in the manuscript "Prolonged delays in human microbiota transmission after a controlled antibiotic perturbation," available here:
https://www.biorxiv.org/content/10.1101/2023.09.26.559480v2

The code here processes raw sequencing reads, identifies shared strains, and performs additional custom analyses.

**System Requirements**

Main software dependencies:
- snakemake version >=5.3
- python version 2.7
- R version >4.0
- MIDAS
- skewer
- bowtie2

The full conda environments used to test and run this pipeline are available here:
https://github.com/ksxue/household-transmission-after-abx/blob/main/envs/skewer-no-builds.yml
https://github.com/ksxue/household-transmission-after-abx/blob/main/envs/bowtie2-no-builds.yml
https://github.com/ksxue/household-transmission-after-abx/blob/main/envs/MIDASpython2-no-builds.yml
https://github.com/ksxue/household-transmission-after-abx/blob/main/envs/Renv-minimal.yml

No non-standard hardware is required.

**Installation Instructions**

The conda environments used to test and run this code are listed above. The typical time to install all dependencies is ~2 hours.

**Demo and Instructions for Use**

To run the main part of the pipeline, execute the following command from the top level of this Github repo:
snakemake

To run the custom analyses, execute the following command from the top level of this Github repo:
workflow/analysis/runAll.R

Expected outputs are included in this Github repo. The expected run time is ~3 hours.




