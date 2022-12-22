#!/usr/bin/env bash

#BSUB -J view[1-3]
#BSUB -o logs/view_%J.out
#BSUB -e logs/view_%J.err
#BSUB -R "select[mem>24] rusage[mem=24] " 
#BSUB -q rna
#BSUB -n 6

module load samtools

SAMPLES=(
HPAP-051
HPAP-053
HPAP-057
)


sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}

GSEA=HPAP

input=./$GSEA/STAR_genome_output/${sample}.bamAligned.out.sam
output=./$GSEA/STAR_genome_output/${sample}.bam
sort_output=./$GSEA/STAR_genome_output/${sample}_sorted.bam

echo $sample

samtools view -bS $input > $output
samtools sort $output -o $sort_output