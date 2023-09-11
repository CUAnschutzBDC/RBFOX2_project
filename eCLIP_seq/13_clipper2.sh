#!/usr/bin/env bash

#BSUB -J clipper[1-8]
#BSUB -o logs/clipper2_%J.out
#BSUB -e logs/clipper2_%J.err
#BSUB -R "select[mem>40] rusage[mem=40] " 
#BSUB -q cranio
#BSUB -n 4


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

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}

INPUT_bam=./output/08_samtools_sort/${sample}_r1.bam
OUTPUT_bed=./output/13_clipper_output/${sample}_r1.peaks.bed

echo $sample

clipper -b $INPUT_bam -o $OUTPUT_bed -s mm10 --processors=4 

# last modified on 7/5/23 to try removing  --binomial=1 --poisson-cutoff=1 for use in replicate analysis

