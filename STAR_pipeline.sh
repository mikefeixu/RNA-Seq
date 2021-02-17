#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -pe smp 8
#$ -N RNA-Seq-STAR
#$ -S /bin/bash
#$ -l h_vmem=8g
#$ -l h_rt=20:00:00
#$ -t 1-9

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
statsdir=$sourcedir/stats
hitcountsdir=$sourcedir/hitcounts
finalcountsdir=$sourcedir/finalcounts
genomedir=/gs/gsfs0/users/xfei/ref/STAR4
gtf=$genomedir/mm10.ERCC92.gtf

mkdir -p $trimdir
mkdir -p $mappingdir
mkdir -p $outtmpdir
mkdir -p $statsdir
mkdir -p $hitcountsdir
mkdir -p $finalcountsdir

# Get sample name from sample_list.txt
sample=$(awk "NR==${SGE_TASK_ID}" $sourcedir/sample_list.txt)
echo "Started task ${SGE_TASK_ID} for ${sample} "; date

# Trim raw data
trim_galore --paired $indir/${sample}_1.fq.gz $indir/${sample}_2.fq.gz -o $trimdir/ --basename ${sample}
echo "Finished trim_galore trimming"; date

#Count reads
zcat $indir/${sample}_1.fq.gz | echo $(($(wc -l)/4)) > $statsdir/${sample}_raw.counts.txt
zcat $trimdir/${sample}_val_1.fq.gz | echo $(($(wc -l)/4)) > $statsdir/${sample}_trimmed.counts.txt
echo "Finished counting of trimmed reads"; date

# Mapping
echo "Start mapping... for sample ${sample}"; date
STAR --runThreadN 8 --genomeDir $genomedir --readFilesIn $trimdir/${sample}_val_1.fq.gz $trimdir/${sample}_val_2.fq.gz --readFilesCommand zcat --outFileNamePrefix $mappingdir/${sample} --outSAMtype BAM Unsorted --outReadsUnmapped Fastx --outTmpDir $outtmpdir/${sample}
samtools sort --threads 8 $mappingdir/${sample}Aligned.out.bam -m 1000000000 -o $bamdir/${sample}.bam
samtools index $bamdir/${sample}.bam

##Check whethr the library is stranded***
#infer_experiment.py -i ${sample}.bam -r mm10_RefSeq.bed

# Feature Counts
#If "gene_id" is preferred, using option "-g gene_id" instead. Same to other types of ids.
featureCounts -a $gtf -o $hitcountsdir/${sample}.counts.txt --largestOverlap -t exon -g gene_name -s 0 -T 8 $bamdir/${sample}.bam
sed -i "s|$bamdir/${sample}.bam|${sample}|g" $hitcountsdir/${sample}.counts.txt
sed -i "s|transcript:||g" $hitcountsdir/${sample}.counts.txt
tail -n +2 $hitcountsdir/${sample}.counts.txt > $finalcountsdir/${sample}.counts.txt
