#!/usr/bin/env bash

#BSUB -J sam_view[1-4]
#BSUB -o logs/sam_view_%J.out
#BSUB -e logs/sam_view_%J.err
#BSUB -R "select[mem>64] rusage[mem=64] " 
#BSUB -q rna
#BSUB -n 6


SAMPLES=(
IgG_CLIP
IgG_INPUT
RBFOX2_CLIP
RBFOX2_INPUT
)

module load samtools

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}

input=./output/08_samtools_sort/${sample}_merged.bam

output=./output/08_samtools_sort/${sample}_r1_merged.bam

echo $sample

samtools view -hb -f 64 $input > $output

# last modified on 7/7/23 for use in merged replicate analysis
