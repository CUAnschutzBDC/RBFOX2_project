#!/usr/bin/env bash

#BSUB -J fastqc1[1-8]
#BSUB -o logs/fastqc1_%J.out
#BSUB -e logs/fastqc1_%J.err
#BSUB -R "select[mem>40] rusage[mem=40] " 
#BSUB -q rna
#BSUB -n 6

module load fastqc/0.11.7

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

fastqc1=./output/01_fastqc

r1_file=./fastq/${sample}_R1_001.fastq.gz 
r2_file=./fastq/${sample}_R2_001.fastq.gz 

fastqc $r1_file -o $fastqc1
fastqc $r2_file -o $fastqc1

# last modified on 6/29/23 for use on the replicate analysis