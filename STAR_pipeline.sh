#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -pe smp 8
#$ -N RNA-Seq-STAR
#$ -S /bin/bash
#$ -l h_vmem=8g
#$ -l h_rt=20:00:00
#$ -t 1-10

# Load modules from HPC
module load cutadapt/1.8.3/python.2.7.8 #Needed by trim_galore
module load trim_galore/0.6.5
module load samtools/1.9/gcc.4.4.7
module load java/1.8.0_20
module load picard
module load STAR
module load subread/1.5.0-p1/gcc.4.4.7
module load RSeQC/2.6.4/python.2.7.8

echo This is task $SGE_TASK_ID

# Make directory structure
sourcedir=$(pwd)
indir=$sourcedir/00_fastq
trimdir=$sourcedir/00_fastq/trimmed
mappingdir=$sourcedir/01_mapping
outtmpdir=$sourcedir/01_mapping/tmp
bamdir=$sourcedir/02_bam
statsdir=$sourcedir/stats
hitcountsdir=$sourcedir/hitcounts
finalcountsdir=$sourcedir/finalcounts
genomedir="Your STAR Index Directory"
gtf=$genomedir/mm10_no_rRNA.gtf
RefSeqbed=$genomedir/mm10_RefSeq_Ensembl.bed

mkdir -p $trimdir
mkdir -p $mappingdir
mkdir -p $outtmpdir
mkdir -p $bamdir
mkdir -p $statsdir
mkdir -p $hitcountsdir
mkdir -p $finalcountsdir

# Get sample name from sample_list.txt
sample=$(awk "NR==${SGE_TASK_ID}" $sourcedir/sample_list.txt)
echo "Started task ${SGE_TASK_ID} for ${sample} "; date

# Trim raw data, for single ended data use commented lines: 
# trim_galore $indir/${sample}.fastq.gz -o $trimdir/ --basename ${sample}
trim_galore --paired $indir/${sample}_1.fq.gz $indir/${sample}_2.fq.gz -o $trimdir/ --basename ${sample}
echo "Finished trim_galore trimming"; date

#Count reads
# zcat $indir/${sample}.fastq.gz | echo $(($(wc -l)/4)) > $statsdir/${sample}_raw.counts.txt
zcat $indir/${sample}_1.fq.gz | echo $(($(wc -l)/4)) > $statsdir/${sample}_raw.counts.txt
# zcat $trimdir/${sample}_trimmed.fq.gz | echo $(($(wc -l)/4)) > $statsdir/${sample}_trimmed.counts.txt
zcat $trimdir/${sample}_val_1.fq.gz | echo $(($(wc -l)/4)) > $statsdir/${sample}_trimmed.counts.txt
echo "Finished counting of trimmed reads"; date

# Mapping
echo "Start mapping... for sample ${sample}"; date
# STAR --runThreadN 8 --genomeDir $genomedir --readFilesIn $trimdir/${sample}_trimmed.fq.gz --readFilesCommand zcat --outFileNamePrefix $mappingdir/${sample} --outSAMtype BAM Unsorted --outReadsUnmapped Fastx --outTmpDir $outtmpdir/${sample}
STAR --runThreadN 8 --genomeDir $genomedir --readFilesIn $trimdir/${sample}_val_1.fq.gz $trimdir/${sample}_val_2.fq.gz --readFilesCommand zcat --outFileNamePrefix $mappingdir/${sample} --outSAMtype BAM Unsorted --outReadsUnmapped Fastx --outTmpDir $outtmpdir/${sample}
samtools sort --threads 8 $mappingdir/${sample}Aligned.out.bam -m 1000000000 -o $bamdir/${sample}.bam
samtools index $bamdir/${sample}.bam

##Check whethr the library is stranded***
infer_experiment.py -i $bamdir/${sample}.bam -r $RefSeqbed
# None Stranded data look like:
# Reading reference gene model mm10_RefSeq.bed ... Done
# Loading SAM/BAM file ...  Total 200000 usable reads were sampled
# This is PairEnd Data
# Fraction of reads failed to determine: 0.0212
# Fraction of reads explained by "1++,1--,2+-,2-+": 0.4974
# Fraction of reads explained by "1+-,1-+,2++,2--": 0.4814

# Feature Counts
#If "gene_id" is preferred, using option "-g gene_id" instead. Same to other types of ids.
featureCounts -a $gtf -o $hitcountsdir/${sample}.counts.txt --largestOverlap -t exon -g gene_name -s 0 -T 8 $bamdir/${sample}.bam
sed -i "s|$bamdir/${sample}.bam|${sample}|g" $hitcountsdir/${sample}.counts.txt
sed -i "s|transcript:||g" $hitcountsdir/${sample}.counts.txt
tail -n +2 $hitcountsdir/${sample}.counts.txt > $finalcountsdir/${sample}.counts.txt
