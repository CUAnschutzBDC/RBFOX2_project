#!/usr/bin/env bash

#BSUB -J fastqc1[1-8]
#BSUB -o logs/fastqc2_%J.out
#BSUB -e logs/fastqc2_%J.err
#BSUB -R "select[mem>40] rusage[mem=40] " 
#BSUB -q rna
#BSUB -n 6

module load fastqc

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

echo $sample

r1_file=./output/03_rmRep/${sample}_rep.bamUnmapped.out.mate1
r2_file=./output/03_rmRep/${sample}_rep.bamUnmapped.out.mate2

fastqc $r1_file -o ./output/04_fastqc2
fastqc $r2_file -o ./output/04_fastqc2


# last modified on 6/30/23 for use in the replicate analysis
