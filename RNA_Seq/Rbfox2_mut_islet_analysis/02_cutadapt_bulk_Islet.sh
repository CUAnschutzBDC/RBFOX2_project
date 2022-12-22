#!/usr/bin/env bash

#BSUB -J fastqc1[1-6]
#BSUB -o logs/fastqc_%J.out
#BSUB -e logs/fastqc_%J.err
#BSUB -R "select[mem>40] rusage[mem=40] " 
#BSUB -q rna
#BSUB -n 6

module load fastqc/0.11.7

SAMPLES=(
822_S34
824_S35
825_S36
826_S37
827_S38
828_S39
)

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}

echo $sample

r1_file=Islet_Rbfox2/fastq/${sample}_L003_R1_001.fastq.gz
r2_file=Islet_Rbfox2/fastq/${sample}_L003_R2_001.fastq.gz

cutadapt -b AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -B AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o fastq/${sample}_L003_R1_001_trimmed.fastq -p fastq/${sample}_L003_R2_001_trimmed.fastq $r1_file $r2_file
