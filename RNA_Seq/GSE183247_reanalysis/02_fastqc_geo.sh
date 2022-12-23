#!/usr/bin/env bash

#BSUB -J fastqc1[1-10]
#BSUB -o logs/fastqc2_%J.out
#BSUB -e logs/fastqc2_%J.err
#BSUB -R "select[mem>40] rusage[mem=40] " 
#BSUB -q rna
#BSUB -n 6

module load fastqc/0.11.7

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

echo $sample

r1_file=${sample}_1.fastq
r2_file=${sample}_2.fastq

fastqc $r1_file -o $fastqc1
fastqc $r2_file -o $fastqc1
