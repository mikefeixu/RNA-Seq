# RNA-Seq User Guide

RNA-Seq Piplines using Tophat or STAR

## Table of Contents

- [RNA-Seq User Guide](#RNA-Seq-user-guide)
  - [Table of Contents](#table-of-contents)
    - [Overview](#running-the-pipeline)
    - [Build STAR Index](#build-star-index)
  - [Report issues or feature requests](#report-issues-or-feature-requests)

### Overview

- Prepare index and data:
  - Build index
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
STAR --runThreadN 8 --runMode genomeGenerate --genomeDir . --genomeFastaFiles mm10.ERCC92.fa --sjdbGTFfile mm10.ERCC92.gtf --sjdbOverhang 149  --genomeChrBinNbits 18 --limitGenomeGenerateRAM 48524399488
```

### Report issues or feature requests

- Open git repository link: [Issues and Feature Requests](https://github.com/mikefeixu/RNA-Seq/issues)
- Click "New Issue"
- Enter the details and submit.