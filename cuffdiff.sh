#!/bin/bash
#$ -cwd
#$ -j y
#$ -N RPKM
#$ -S /bin/bash
#$ -pe smp 6
#$ -cwd
#$ -l h_vmem=20g
#$ -hold_jid 5514137

module load cufflinks/2.2.1/gcc.4.4.7
echo "Started "; date

sourcedir=$(pwd)
bamdir=$sourcedir/01_mapping
genomedir=/gs/gsfs0/users/xfei/ref/Bowtie2/mm10
gtf=mm10_refSeq_Mar132017_repIDwithName_ERCC.gtf

mkdir -p $bamdir
cd $bamdir

cuffdiff $genomedir/$gtf 1_Rinf_WT_1_day0_accepted_hits.bam 2_Rinf_WT_2_day0_accepted_hits.bam 3_Rinf_KO_1_day0_accepted_hits.bam 4_Rinf_KO_2_day0_accepted_hits.bam 5_Rinf_WT_1_day3_accepted_hits.bam 6_Rinf_WT_2_day3_accepted_hits.bam 7_Rinf_KO_1_day3_accepted_hits.bam 8_Rinf_KO_2_day3_accepted_hits.bam 9_Rinf_WT_1_day6_accepted_hits.bam 10_Rinf_WT_2_day6_accepted_hits.bam 11_Rinf_KO_1_day6_accepted_hits.bam 12_Rinf_KO_2_day6_accepted_hits.bam
#cuffdiff $genomedir/$gtf 1_Rinf_WT_1_day0_accepted_hits.bam,2_Rinf_WT_2_day0_accepted_hits.bam 3_Rinf_KO_1_day0_accepted_hits.bam,4_Rinf_KO_2_day0_accepted_hits.bam 5_Rinf_WT_1_day3_accepted_hits.bam,6_Rinf_WT_2_day3_accepted_hits.bam 7_Rinf_KO_1_day3_accepted_hits.bam,8_Rinf_KO_2_day3_accepted_hits.bam 9_Rinf_WT_1_day6_accepted_hits.bam,10_Rinf_WT_2_day6_accepted_hits.bam 11_Rinf_KO_1_day6_accepted_hits.bam,12_Rinf_KO_2_day6_accepted_hits.bam

echo "Finished "; date
