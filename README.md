# RNA-Seq User Guide

RNA-Seq Piplines using Tophat or STAR

## Table of Contents

- [RNA-Seq User Guide](#RNA-Seq-user-guide)
  - [Table of Contents](#table-of-contents)
    - [Overview](#running-the-pipeline)
    - [Build STAR Index](#build-star-index-locally)
    - [Build STAR Index on HPC](#build-star-index-on-HPC)
    - [Download data](#download-data)
    - [Run STAR pipeline](#run-star-pipeline)
    - [Run fastqc and RSeQC](#run-fastqc-and-RSeQC)
    - [Run DEG Analysis](#run-deg-analysis)
  - [Report issues or feature requests](#report-issues-or-feature-requests)

### Overview

- Prepare index and data:
  - Build STAR index
  - Download data
- STAR Pipeline:
  - Run STAR_pipeline.sh on HPC environment
  - Run deseq2STAR.R
- Data Visualization:
  - Run BubblePlotGoTerm.R
  - Run HeatMap.R
  - Run VennDiagram.R
  - Run BarPlot.R

#### Build STAR index locally

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

#### Build STAR index on HPC

- Submit STAR_index.sh

```bash
# Config the genomedir
# Update the genome and gtf file path for download
qsub STAR_index.sh

# Check the download logs to ensure no transfer errors. 
```

#### Download data

- Download published data with HPC environment. Download scripts and example are available at [SRA Download](https://github.com/mikefeixu/sra-download)

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
# Add sample_list.txt to script directory and update it with sample names

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

#### Run DEG analysis

- Run deseq2STAR.R

```r
# Config comparison and projectdir before running the script.
# Create a folder for each comparion, add the counts files (output of STAR_pipeline.sh) for samples from both groups to this folder
# Add Conditions.txt to each comparion folder, and update sample, group, and batch info for according comparions. Leave batch info blank if no batch effect expected.
# Add SelectedLabledGenes.csv to Selected_genes folder, and update it with the genes you want to lable in the MA plot and Volcano plot
# Install below packages if not installed
# Recommend to run it in Rstudio
# Run deseq2STAR.R for each comparison first, then run the R scripts in plots folder for BarPlot, BubblePlot, Heatmap, and VennDiagrams.
# To included all samples in the PCA plot, run a pseduo comparison that includes all of the samples.

library(DESeq2)
library(gtools)
library(dplyr)
library(RColorBrewer)
library(pheatmap)
library(ggplot2)

# Run below scripts in Rstudio
# deseq2STAR.R
# VennDiagram.R
# HeatMap.R
# BarPlot.R
# BubblePlotGOTerm.R
```

### Report issues or feature requests

- Open git repository link: [Issues and Feature](https://github.com/mikefeixu/RNA-Seq/issues)
- Click "New Issue"
- Enter the details and submit.
