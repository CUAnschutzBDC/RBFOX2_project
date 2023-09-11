#!/usr/bin/env bash

#BSUB -J sam_sort[1-8]
#BSUB -o logs/sam_sort_%J.out
#BSUB -e logs/sam_sort_%J.err
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

sort_dir=./output/08_samtools_sort

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}

input=./output/07_barcode_collapse/${sample}_collapse.bam

output=$sort_dir/${sample}_collapse_sorted.bam

echo $sample

samtools sort $input -o $output

# last modified 6/30/23 for use in the replicate analysis
