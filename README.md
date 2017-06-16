# rnaseq-pipeline

## Step1: Alignment
STAR --genomeDir [path-to-genome-index] --readFilesIn Sample_R1.fastq.gz Sample_R2.fastq.gz --readFilesCommand zcat --runThreadN 4
# This will output bam files. 

## Step2: quantification
featureCounts -a genes.saf -o output.counts [BAM FILES] -F SAF -T 8


## quality checks
fastqc *.fastq.gz # check fastq files
samtools flagstat [BAM FILES] # check 
