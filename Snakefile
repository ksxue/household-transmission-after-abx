from os.path import join
import pandas as pd
import os

configfile: "config/config.yaml"

# Convert list of samples to a dataframe.
df=pd.read_table(config['sampleFileRaw'], header=None)
df.columns=["sample","samplelane","read1","read2","trim1","trim2"]

# Parse subject and household from sample names.
df['samplename']=df['sample']
df[['subject','timepoint']]=df.samplename.str.split("-",expand=True)
df['household']=df['subject'].astype(str).str[0:2]
df['studyArm']=df['subject'].astype(str).str[0:1]

# Generate lists of subjects and households.
subjects=list(set(df['subject'].tolist()))
households=list(set(df['household'].tolist()))
samplelanes=list(set(df['samplelane'].tolist()))
samples=list(set(df['samplename'].tolist()))

# Parse the list of species analyzed in the MIDAS snps module.
if(os.path.isdir("workflow/out/midasOutput/snps")):
    # Iterate through the species analyzed by the snps module.
    dirs=[x[0] for x in os.walk("workflow/out/midasOutput/snps/HouseholdTransmission-Stool")]
    # Remove the first element, which is the directory itself without subdirectories.
    dirs.pop(0)
    # Parse the species names and generate a list.
    snpsSpecies=[]
    for species in dirs:
        snpsSpecies.append(species.split("/")[5])
        
# Import the list of species of special interest in household XB.
dfXB=pd.read_table(config['XBspeciesOfInterest'], header=None)
dfXB.columns=["species"]
XBspeciesOfInterest=list(set(dfXB['species'].tolist()))

rule all:
    input:
        #expand("/scratch/groups/relman/kxue/household-transmission-mgx/trimmed/HouseholdTransmission-Stool-{sample}-trimmed-pair1.fastq.gz",sample=samplelanes),
        #expand("/scratch/groups/relman/kxue/household-transmission-mgx/trimmed/HouseholdTransmission-Stool-{sample}.done",sample=samplelanes)
        #expand("workflow/out/filter/HouseholdTransmission-Stool-{sample}-filtered.1.fastq.gz",sample=samples),
        #expand("workflow/out/midasOutput/HouseholdTransmission-Stool-{sample}/species/species_profile.txt",sample=samples),
        #"workflow/out/midasOutput/species/species_profile_all.txt"
        #expand("workflow/out/midasOutput/species/abundantSpecies_{subject}.txt", subject=subjects),
        #expand("workflow/out/midasOutput/species/abundantSpecies_{household}.txt", household=households)
        #expand("workflow/out/midasOutput/HouseholdTransmission-Stool-{sample}/snps/summary.txt", sample=samples)
        #expand("workflow/out/midasOutput/snps/HouseholdTransmission-Stool/{species}/snps_freq.txt.gz", species=snpsSpecies)
        #expand("workflow/report/calculateFixedDifferences/{species}/done.txt", species=snpsSpecies),
        #"workflow/report/calculateFixedDifferences/Bacteroides_stercoris_56735/done.txt"
        #expand("workflow/report/performStrainFishing/{species}/done.txt", species=snpsSpecies)
        #expand("workflow/report/calculateFixedDifferences/{species}/fixedDiffs_filterCoverage.txt.gz", species=snpsSpecies)
        #"workflow/report/performStrainFishing/Bacteroides_stercoris_56735/done.txt"
        #"workflow/report/performStrainFishing/Akkermansia_muciniphila_55290/done.txt"
        #"workflow/report/performStrainFishing/removeEmptyStrainFishing.done"
        #"workflow/out/T6SS/HouseholdTransmission-Stool-XBA-001.bam",
        #"workflow/out/T6SS/HouseholdTransmission-Stool-XBA-058.bam",
        #"workflow/out/T6SS/HouseholdTransmission-Stool-XBA-519.bam"
        #expand("workflow/out/T6SS/HouseholdTransmission-Stool-{sample}.bam", sample=samples),
        #"workflow/out/T6SS/T6SS-coverage.txt",
        #join(config["T6SSdir"],"T6SS-totalReads.txt")
        expand("workflow/out/midasOutput/snps/HouseholdTransmission-Stool/{species}/XB.rds", species=XBspeciesOfInterest)
        


#include: "workflow/rules/processRawReads.smk"
#include: "workflow/rules/runMIDAS.smk"
include: "workflow/rules/postMIDAS.smk"
include: "workflow/rules/mapReadsT6SS.smk"
