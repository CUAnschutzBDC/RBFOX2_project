#!/usr/bin/env bash

#BSUB -J sam_index[1-8]
#BSUB -o logs/sam_index_%J.out
#BSUB -e logs/sam_index_%J.err
#BSUB -R "select[mem>64] rusage[mem=64] " 
#BSUB -q rna
#BSUB -n 6


SAMPLES=(
IgG_CLIP_S25_L001
IgG_IN_S27_L001
Rbfox2_CLIP_S24_L001
Rbfox2_IN_S26_L001
IgG_CLIP_S30_L002
IgG_INPUT_S28_L002
RBFOX2_CLIP_S31_L002
RBFOX2_INPUT_S29_L002
)

module load samtools

index_dir=./output/08_samtools_sort

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}

input=./output/08_samtools_sort/${sample}_collapse_sorted.bam
output=$index_dir/${sample}.bam.bai

echo $sample

samtools index $input $output

# last modified on 6/30/23 for use in the replicate analysis
