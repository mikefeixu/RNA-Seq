#!/bin/bash
#
# task.sh
#
#$ -cwd
#$ -j y
#$ -pe smp 8
#$ -N RNA-Seq
#$ -S /bin/bash
#$ -l h_vmem=8g
#$ -l h_rt=20:00:00
#$ -t 1-12
#$ -M fei.xu@einsteinmed.org

module load STAR
module load samtools
module load subread/1.5.0-p1/gcc.4.4.7
echo This is task $SGE_TASK_ID

sourcedir=$(pwd)
mappingdir=$sourcedir/01_mapping
indir=$sourcedir/00_fastq/trimmed
outtmpdir=$sourcedir/01_mapping/tmp
bamdir=$sourcedir/bam
hitcountsdir=$sourcedir/hitcounts
finalcountsdir=$sourcedir/finalcounts
mkdir -p $mappingdir
mkdir -p $outtmpdir
mkdir -p $bamdir
mkdir -p $hitcountsdir
mkdir -p $finalcountsdir
#cd $hitcountsdir

sample=$(awk "NR==${SGE_TASK_ID}" $sourcedir/sample_list.txt)

echo $sample
STAR --runThreadN 8 --genomeDir /gs/gsfs0/users/xfei/ref/STAR3 --readFilesIn $indir/${sample}.fastq.gz --readFilesCommand zcat --outFileNamePrefix $mappingdir/${sample} --outSAMtype BAM Unsorted --outReadsUnmapped Fastx --outTmpDir $outtmpdir/${sample}
samtools sort --threads 8 $mappingdir/${sample}Aligned.out.bam -m 1000000000 -o $bamdir/${sample}.bam
samtools index $bamdir/${sample}.bam

##Check whethr the library is stranded***
featureCounts -a /gs/gsfs0/users/xfei/ref/STAR2/mm10.ERCC92.gtf -o $hitcountsdir/${sample}.counts.txt --largestOverlap -t exon -g gene_name -s 0 -T 8 $bamdir/${sample}.bam
sed -i "s|$bamdir/${sample}.bam|${sample}|g" $hitcountsdir/${sample}.counts.txt
sed -i "s|transcript:||g" $hitcountsdir/${sample}.counts.txt
tail -n +2 $hitcountsdir/${sample}.counts.txt > $finalcountsdir/${sample}.counts.txt