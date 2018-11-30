# rnaseq-pipeline
## File requirement.
First make your project directory, and make a fastq directory in the project folder. `project/fastq`. 
Put all of your fastq files in the `project/fastq` folder. Fastq files could be in .fastq or .gz or .bz2 format. 
Name your fastq files either as *.fastq.[gz/bz2] (single read) or *_R1.fastq and *_R2.fastq (paired-end reads). 

### Step1: Alignment
STAR --genomeDir [path-to-genome-index] --readFilesIn Sample_R1.fastq.gz Sample_R2.fastq.gz --readFilesCommand zcat --runThreadN 4
# This will output bam files. 

### Step2: quantification
featureCounts -a genes.saf -o output.counts [BAM FILES] -F SAF -T 8


### quality checks
fastqc *.fastq.gz # check fastq files
samtools flagstat [BAM FILES] # check 
