#!/bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 8
#$ -N RNA-Seq-Tophat
#$ -S /bin/bash
#$ -l h_vmem=16g
#$ -l h_rt=20:00:00
#$ -t 1-9


module load cutadapt/1.8.3/python.2.7.8       # Needed by trim_galore
module load trim_galore/0.6.5
module load java/1.8.0_20
module load tophat
module load HTSeq/0.6.1/python.2.7.8
module load samtools/1.9/gcc.4.4.7

echo This is task $SGE_TASK_ID

sourcedir=$(pwd)
indir=$sourcedir/00_fastq
trimdir=$sourcedir/00_fastq/trimmed
statsdir=$sourcedir/stats
mappingdir=$sourcedir/01_mapping
htseqresults=$sourcedir/htseqresults
genomedir=/gs/gsfs0/users/xfei/ref/Bowtie2/mm10
gtf=mm10_refSeq_Mar132017_repIDwithName_ERCC.gtf

echo This is task $SGE_TASK_ID

mkdir -p $trimdir
mkdir -p $statsdir
mkdir -p $mappingdir
mkdir -p $htseqresults

sample=$(awk "NR==${SGE_TASK_ID}" $sourcedir/sample_list.txt)
echo "Started task ${SGE_TASK_ID} for sample: ${sample}"; date

# Trim raw data
trim_galore --paired $indir/${sample}_1.clean.fq.gz $indir/${sample}_2.clean.fq.gz -o $trimdir/

#Count reads
zcat $indir/${sample}_1.clean.fq.gz | echo $((`wc -l`/4)) > $statsdir/${sample}_raw.counts.txt
zcat $trimdir/${sample}_1.clean_val_1.fq.gz | echo $((`wc -l`/4)) > $statsdir/${sample}_trimmed.counts.txt

echo "mapping..."; date
# Map to genome
tophat -p 4 -G $genomedir/$gtf -o $mappingdir/${sample}_tophatout $genomedir/mm10 $trimdir/${sample}_1.clean_val_1.fq.gz $trimdir/${sample}_2.clean_val_2.fq.gz 

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
