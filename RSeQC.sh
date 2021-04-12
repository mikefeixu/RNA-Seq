#!/bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 4
#$ -N RSeQC
#$ -S /bin/bash
#$ -l h_vmem=16g
#$ -l h_rt=20:00:00
#$ -t 1-10

module load RSeQC/2.6.4/python.2.7.8
echo This is task $SGE_TASK_ID; date

sourcedir=$(pwd)
trimmdir=$sourcedir/00_fastq/trimmed
bamdir=$sourcedir/02_bam
clipdir=$sourcedir/02_bam/clip
geneBodydir=$sourcedir/02_bam/geneBody
junctiondir=$sourcedir/02_bam/junction
RPKMdir=$sourcedir/02_bam/RPKM
genomedir=/gs/gsfs0/users/xfei/ref/STAR
RefSeqbed=$genomedir/mm10_RefSeq_Ensembl.bed

mkdir -p $clipdir
mkdir -p $geneBodydir
mkdir -p $junctiondir
mkdir -p $RPKMdir

sample=$(awk "NR==${SGE_TASK_ID}" $sourcedir/sample_list.txt)

# FastQC
echo "Started task ${SGE_TASK_ID} for sample: ${sample} FastQC";
cd $trimmdir
fastqc $trimmdir/${sample}_val_1.fq.gz $trimmdir/${sample}_val_2.fq.gz;
echo "Finished sample: ${sample} FastQC";

# RSeQC
cd $bamdir
sample=$(awk "NR==${SGE_TASK_ID}" $sourcedir/sample_list.txt)
echo "Started task ${SGE_TASK_ID} for sample: ${sample} RSeQC";
clipping_profile.py -i ${sample}.bam -o $clipdir/${sample} -s PE
geneBody_coverage.py -i ${sample}.bam -o $geneBodydir/${sample} -r $RefSeqbed -l 2000
junction_annotation.py -i ${sample}.bam -o $junctiondir/${sample} -r $RefSeqbed
RPKM_saturation.py -i ${sample}.bam -o $RPKMdir/${sample} -r $RefSeqbed
echo "Finished sample: ${sample} RSeQC";
