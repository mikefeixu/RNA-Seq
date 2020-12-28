#!/bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 8
#$ -N RNA-Seq-R2
#$ -S /bin/bash
#$ -l h_vmem=16g
#$ -l h_rt=20:00:00
#$ -t 1-12

module load java/1.8.0_20
module load tophat
module load HTSeq/0.6.1/python.2.7.8
module load samtools/1.9/gcc.4.4.7

echo This is task $SGE_TASK_ID

sourcedir=$(pwd)
mappingdir=$sourcedir/01_mapping
htseqresults=$sourcedir/htseqresults
indir=/gs/gsfs0/users/xfei/Qin/RNA-Seq/00_fastq/trimmed
genomedir=/gs/gsfs0/users/xfei/ref/Bowtie2/mm10
gtf=mm10_refSeq_Mar132017_repIDwithName_ERCC.gtf

mkdir -p $mappingdir
mkdir -p $htseqresults

sample=$(awk "NR==${SGE_TASK_ID}" $sourcedir/sample_list.txt)
echo "Started task ${SGE_TASK_ID} for sample: ${sample} Renaming"; date

# Mapping
tophat -p 4 -G $genomedir/$gtf -o $mappingdir/${sample}_tophatout $genomedir/mm10 $indir/${sample}.fastq.gz

# Index the bam
samtools index $mappingdir/${sample}_tophatout/accepted_hits.bam

# Collect the stats
samtools flagstat $mappingdir/${sample}_tophatout/accepted_hits.bam > $mappingdir/${sample}_tophatout/accepted_hits.flastats

# HTseq Count
htseq-count -f bam --stranded=no -r name $mappingdir/${sample}_tophatout/accepted_hits.bam $genomedir/$gtf > $htseqresults/${sample}.htseq.txt

# Convert to sam file
samtools view -h $mappingdir/${sample}_tophatout/accepted_hits.bam > $mappingdir/${sample}_accepted_hits.sam

# Rename
mv $mappingdir/${sample}_tophatout/accepted_hits.bam $mappingdir/${sample}_accepted_hits.bam
mv $mappingdir/${sample}_tophatout/accepted_hits.bam.bai $mappingdir/${sample}_accepted_hits.bam.bai

echo "Finished Renaming"; date
