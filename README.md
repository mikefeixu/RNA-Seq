# RNA-Seq User Guide

RNA-Seq Piplines using Tophat or STAR

## Table of Contents

- [RNA-Seq User Guide](#RNA-Seq-user-guide)
  - [Table of Contents](#table-of-contents)
    - [Overview](#running-the-pipeline)
    - [Build STAR Index](#build-star-index)
    - [Download data](#download-data)
    - [Run STAR pipeline](#run-star-pipeline)
  - [Report issues or feature requests](#report-issues-or-feature-requests)

### Overview

- Prepare index and data:
  - Build STAR index
  - Download data
- STAR Pipeline:
  - Run STAR_pipeline.sh on HPC environment
  - Run deseq2STAR.R
- Tophat Pipeline:
  - Run Tophat_pipeline.sh on HPC environment
  - Run deseq2Tophat.R
- Data Visualization:
  - Run BubblePlotGoTerm.R
  - Run HeatMap.R
  - Run VennDiagram.R

#### Build STAR index

- Download genome

```bash
wget http://ftp.ensembl.org/pub/release-103/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
```

- Download annotation file

```bash
wget http://ftp.ensembl.org/pub/release-103/gtf/mus_musculus/Mus_musculus.GRCm39.103.gtf.gz
```

- Preprocess genome files

```bash
gunzip *.gz
mv Mus_musculus.GRCm39.dna.primary_assembly.fa mm10.fa
mv Mus_musculus.GRCm39.103.gtf mm10.gtf
grep -v -i 'gene_biotype "rRNA";\|gene_biotype "Mt_rRNA";' mm10.gtf > mm10_no_rRNA.gtf
```

- Index genome

```bash
samtools faidx mm10.fa
java -jar -Xmx16g /public/apps/picard/2.17.1/picard.jar CreateSequenceDictionary R=mm10.fa O=mm10.dict
```

- Create STAR Index

```bash
# If read length is 150, use 149 (150-1) for parameter --sjdbOverhang
STAR --runThreadN 8 --runMode genomeGenerate --genomeDir . --genomeFastaFiles mm10.ERCC92.fa --sjdbGTFfile mm10.ERCC92.gtf --sjdbOverhang 149  --genomeChrBinNbits 18 --limitGenomeGenerateRAM 48524399488
```

#### Download data

- Download published data with HPC environment.

```bash
# Config HPC settings in sra-download.sh. Make sure the task number equals to the sample number. (eg. -t 1~9 means 9 sampels for download)
# Update the rename step accordingly for single ended data. Then submit the array jobs to HPC.

qsub sra-download.sh

# Check the download logs to ensure no transfer errors. Move downloaded data to indir in STAR_pipeline.sh after download jobs are completed.
```

#### Run STAR pipeline

- Run STAR pipleine with HPC environment

```bash
# Config HPC settings in sra-download.sh. Make sure the task number equals to the sample number. eg. '-t 1~12' means 12 sampels for download; '-pe smp 8' means using 8 cores of parallel environment(pe) "smp". You need to change it to your pe.

# Configure the directories in the job script. Especially for indir, genomedir, gtf, and RefSeqbed (RefSeqbed is used by infer_experiment.py(RSeQC) for checking library strand info)
# sourcedir=$(pwd)
# indir=$sourcedir/00_fastq
# genomedir=<STAR4 directory>
# gtf=$genomedir/mm10.gtf
# RefSeqbed=$genomedir/mm10_RefSeq_Ensembl.bed

qsub STAR_pipeline.sh

# Check whethr the library is stranded and set the '-s 1'. Default is '-s 0' for non-stranded library.
# infer_experiment.py -i ${sample}.bam -r mm10_RefSeq.bed
# None Stranded data look like:
# Reading reference gene model mm10_RefSeq.bed ... Done
# Loading SAM/BAM file ...  Total 200000 usable reads were sampled
# This is PairEnd Data
# Fraction of reads failed to determine: 0.0212
# Fraction of reads explained by "1++,1--,2+-,2-+": 0.4974
# Fraction of reads explained by "1+-,1-+,2++,2--": 0.4814

```

#### Run fastqc and RSeQC

- Run fastqc and RSeQC pipleine with HPC environment

```bash
# Config HPC settings in sra-download.sh. Make sure the task number equals to the sample number. eg. '-t 1~12' means 12 sampels for download; '-pe smp 8' means using 8 cores of parallel environment(pe) "smp". You need to change it to your pe.

# Configure the directories in the job script. Especially for trimmdir, bamdir (generated or copied from STAR_pipeline.sh), genomedir, and RefSeqbed
# sourcedir=$(pwd)
# trimmdir=$sourcedir/00_fastq/trimmed
# bamdir=$sourcedir/02_bam
# genomedir=<STAR4 directory>
# gtf=$genomedir/mm10.gtf
# RefSeqbed=$genomedir/mm10_RefSeq_Ensembl.bed

qsub RseQC.sh

# Check fastqc and RseQC reports and figures in according folders.
```

### Report issues or feature requests

- Open git repository link: [Issues and Feature](https://github.com/mikefeixu/RNA-Seq/issues)
- Click "New Issue"
- Enter the details and submit.
