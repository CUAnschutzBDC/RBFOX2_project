#!/usr/bin/env bash

#BSUB -J align[1-8]
#BSUB -o logs/align_bulk_%J.out
#BSUB -e logs/align_bulk_%J.err
#BSUB -R "select[mem>24] rusage[mem=24] " 
#BSUB -q rna
#BSUB -n 6

module load STAR/2.7.9a

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

STAR_dir=MIN6_RBFOX2_KD/STAR_genome_output

STAR_file=$STAR_dir/${sample}.bam

echo $sample

r1_file=MIN6_RBFOX2_KD/fastq/${sample}_L002_R1_001_trimmed.fastq
r2_file=MIN6_RBFOX2_KD/fastq/${sample}_L002_R2_001_trimmed.fastq

STAR --genomeDir ../GEO_Datasets/ref/star/mouse/Grm38/index --readFilesIn $r1_file $r2_file --outFileNamePrefix $STAR_file
