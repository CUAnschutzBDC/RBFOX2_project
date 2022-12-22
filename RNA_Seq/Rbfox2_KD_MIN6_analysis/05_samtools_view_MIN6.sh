#!/usr/bin/env bash

#BSUB -J view[1-8]
#BSUB -o logs/view_%J.out
#BSUB -e logs/view_%J.err
#BSUB -R "select[mem>24] rusage[mem=24] " 
#BSUB -q rna
#BSUB -n 6

module load samtools

SAMPLES=(
NT1_S5
NT3_S6
NT5_S7
NT6_S8
RbFox2_1_S1
RbFox2_3_S2
RbFox2_5_S3
RbFox2_6_S4
)

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}

input=MIN6_RBFOX2_KD/STAR_genome_output/${sample}.bamAligned.out.sam
output=MIN6_RBFOX2_KD/${sample}.bam
sort_output=MIN6_RBFOX2_KD/${sample}_sorted.bam



echo $sample

samtools view -bS $input > $output

samtools sort $output -o $sort_output
