#!/usr/bin/env bash

#BSUB -J align[1-6]
#BSUB -o logs/align_bulk_%J.out
#BSUB -e logs/align_bulk_%J.err
#BSUB -R "select[mem>24] rusage[mem=24] " 
#BSUB -q rna
#BSUB -n 6

module load STAR/2.7.9a

SAMPLES=(
822_S34
824_S35
825_S36
826_S37
827_S38
828_S39
)

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}

STAR_dir=Islet_Rbfox2/STAR_genome_output

STAR_file=$STAR_dir/${sample}.bam

echo $sample

r1_file=Islet_Rbfox2/fastq/${sample}_L003_R1_001_trimmed.fastq
r2_file=Islet_Rbfox2/fastq/${sample}_L003_R2_001_trimmed.fastq

STAR --genomeDir ../GEO_Datasets/ref/star/mouse/Grm38/index --readFilesIn $r1_file $r2_file --outFileNamePrefix $STAR_file
