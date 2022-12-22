#!/usr/bin/env bash

#BSUB -J view[1-6]
#BSUB -o logs/view_%J.out
#BSUB -e logs/view_%J.err
#BSUB -R "select[mem>24] rusage[mem=24] " 
#BSUB -q rna
#BSUB -n 6

module load samtools

SAMPLES=(
822_S34
824_S35
825_S36
826_S37
827_S38
828_S39
)

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}

input=./STAR_genome_output/${sample}.bamAligned.out.sam
output=${sample}.bam
sort_output=${sample}_sorted.bam



echo $sample

samtools view -bS $input > $output

samtools sort $output -o $sort_output
