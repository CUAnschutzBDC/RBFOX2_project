#!/usr/bin/env bash

#BSUB -J sort[1-10]
#BSUB -o logs/sort_%J.out
#BSUB -e logs/sort_%J.err
#BSUB -R "select[mem>24] rusage[mem=24] " 
#BSUB -q rna
#BSUB -n 6

module load samtools

SAMPLES=(
SRR15697732
SRR15697733
SRR15697734
SRR15697735
SRR15697736
SRR15697737
SRR15697738
SRR15697739
SRR15697740
SRR15697741
)

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}

input=${sample}.bam
output=${sample}_sorted.bam


echo $sample

samtools sort $input -o $output

sorted_sam=${sample}_sorted.sam

samtools view -h -o $sorted_sam $output