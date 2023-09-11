#!/usr/bin/env bash

#BSUB -J clipper[1-4]
#BSUB -o logs/clipper2_%J.out
#BSUB -e logs/clipper2_%J.err
#BSUB -R "select[mem>40] rusage[mem=40] " 
#BSUB -q cranio
#BSUB -n 4


SAMPLES=(
IgG_CLIP
IgG_INPUT
RBFOX2_CLIP
RBFOX2_INPUT
)

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}

INPUT_bam=./output/08_samtools_sort/${sample}_r1_merged.bam
OUTPUT_bed=./output/13_clipper_output/${sample}_r1_merged.peaks.bed

echo $sample

clipper -b $INPUT_bam -o $OUTPUT_bed -s mm10 --processors=4 

# last modified on 7/7/23 for use in merged replicate analysis

