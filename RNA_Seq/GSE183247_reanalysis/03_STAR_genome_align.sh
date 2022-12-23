#!/usr/bin/env bash

#BSUB -J align[1-10]
#BSUB -o logs/align_%J.out
#BSUB -e logs/align_%J.err
#BSUB -R "select[mem>24] rusage[mem=24] " 
#BSUB -q rna
#BSUB -n 6

module load STAR/2.7.9a

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

STAR_dir=STAR_genome_output

STAR_file=$STAR_dir/${sample}.bam

echo $sample

r1_file=${sample}_1.fastq
r2_file=${sample}_2.fastq

STAR --genomeDir ./ref/star/mouse/Grm38/index --readFilesIn $r1_file $r2_file --outFileNamePrefix $STAR_file
