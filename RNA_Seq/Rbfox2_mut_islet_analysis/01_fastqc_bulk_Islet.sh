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

fastqc $r1_file -o $fastqc1
fastqc $r2_file -o $fastqc1
