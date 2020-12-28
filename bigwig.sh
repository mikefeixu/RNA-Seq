#!/bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 1
#$ -N ConvertToBW
#$ -S /bin/bash
#$ -l h_vmem=16g
#$ -l h_rt=20:00:00
#$ -t 1-12

module load deeptools

echo This is task $SGE_TASK_ID

sourcedir=$(pwd)
mappingdir=$sourcedir/01_mapping

sample=$(awk "NR==${SGE_TASK_ID}" $sourcedir/sample_list.txt)
echo "Started task ${SGE_TASK_ID} for sample: ${sample} convertig to bigwig"; date

# Convert bam to bigwig file
bamCoverage -b $mappingdir/${sample}_accepted_hits.bam -o ${sample}.bw

echo "Finished converting ${sample}.bam to bigwig file"; date
