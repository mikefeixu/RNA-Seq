#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -pe smp 8
#$ -N STAR_index
#$ -S /bin/bash
#$ -l h_vmem=8g
#$ -l h_rt=20:00:00
#$ -t 1

module load samtools
module load picard
module load HTSeq
module load STAR

genomedir="<Your STAR index directory>"
gtf=$genomedir/mm10_no_rRNA.gtf

cd $genomedir||exit

# Update below address to your genome files found on Ensembl, UCSC, or elsewhere
wget http://ftp.ensembl.org/pub/release-103/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
wget http://ftp.ensembl.org/pub/release-103/gtf/mus_musculus/Mus_musculus.GRCm39.103.gtf.gz

# Preprocess genome files
gunzip *.gz
mv Mus_musculus.GRCm39.dna.primary_assembly.fa mm10.fa
mv Mus_musculus.GRCm39.103.gtf mm10.gtf
grep -v -i 'gene_biotype "rRNA";\|gene_biotype "Mt_rRNA";' mm10.gtf > mm10_no_rRNA.gtf

# Index genome
samtools faidx mm10.fa
java -jar -Xmx16g /public/apps/picard/2.17.1/picard.jar CreateSequenceDictionary R=mm10.fa O=mm10.dict

# Create STAR Index
# If read length is 150, use 149 (150-1) for parameter --sjdbOverhang
STAR --runThreadN 8 --runMode genomeGenerate --genomeDir $genomedir --genomeFastaFiles mm10.fa --sjdbGTFfile $gtf --sjdbOverhang 149  --genomeChrBinNbits 18 --limitGenomeGenerateRAM 48524399488
