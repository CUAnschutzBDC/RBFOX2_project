#!/usr/bin/env bash

#BSUB -J align[1-1]
#BSUB -o logs/align_%J.out
#BSUB -e logs/align_%J.err
#BSUB -R "select[mem>24] rusage[mem=24] " 
#BSUB -q rna
#BSUB -n 6

module load STAR/2.7.9a

SAMPLES=(
HPAP-051
)

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}

STAR_dir=HPAP/STAR_genome_output

STAR_file=$STAR_dir/${sample}.bam

echo $sample

r1_file=HPAP/${sample}_mRNASeq_beta_R1_fastq-data.fastq.gz
r2_file=HPAP/${sample}_mRNASeq_beta_R2_fastq-data.fastq.gz

STAR --genomeDir ./ref/star/human/GRCh38/index --readFilesIn $r1_file $r2_file --readFilesCommand zcat --outFileNamePrefix $STAR_file
