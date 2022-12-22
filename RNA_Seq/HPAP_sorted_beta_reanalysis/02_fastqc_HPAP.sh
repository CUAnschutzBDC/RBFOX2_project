#!/usr/bin/env bash

#BSUB -J fastqc1[1-3]
#BSUB -o logs/fastqc2_%J.out
#BSUB -e logs/fastqc2_%J.err
#BSUB -R "select[mem>40] rusage[mem=40] " 
#BSUB -q rna
#BSUB -n 6

module load fastqc/0.11.7

SAMPLES=(
HPAP-051
HPAP-053
HPAP-057
)

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}

echo $sample

fastqc1=HPAP/fastqc
r1_file=HPAP/${sample}_mRNASeq_beta_R1_fastq-data.fastq.gz
r2_file=HPAP/${sample}_mRNASeq_beta_R2_fastq-data.fastq.gz

fastqc $r1_file -o $fastqc1
fastqc $r2_file -o $fastqc1
