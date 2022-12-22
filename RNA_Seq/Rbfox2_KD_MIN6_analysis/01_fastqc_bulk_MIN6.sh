#!/usr/bin/env bash

#BSUB -J fastqc1[1-8]
#BSUB -o logs/fastqc_MIN6%J.out
#BSUB -e logs/fastqc_MIN6%J.err
#BSUB -R "select[mem>40] rusage[mem=40] " 
#BSUB -q rna
#BSUB -n 6

module load fastqc/0.11.7

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

echo $sample

fastqc1=MIN6_RBFOX2_KD/fastqc
r1_file=MIN6_RBFOX2_KD/fastq/${sample}_L002_R1_001.fastq.gz
r2_file=MIN6_RBFOX2_KD/fastq/${sample}_L002_R2_001.fastq.gz

fastqc $r1_file -o $fastqc1
fastqc $r2_file -o $fastqc1
