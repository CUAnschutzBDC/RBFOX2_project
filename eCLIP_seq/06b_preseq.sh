#!/usr/bin/env bash

#BSUB -J preseq[1-8]
#BSUB -o logs/06b_preseq_%J.out
#BSUB -e logs/06b_preseq_%J.err
#BSUB -R "select[mem>40] rusage[mem=40] " 
#BSUB -q rna
#BSUB -n 6

module load samtools
module load bedtools


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


input=./output/06_STAR_genome_output/${sample}.bamAligned.out.bam
sorted_bam=./output/06_STAR_genome_output${sample}_sorted.bam
sorted_bed=./output/06b_preseq/${sample}_sorted.bed
preseq=./output/06b_preseq/${sample}_complecity_output.txt

echo $sample

samtools sort $input -o $sorted_bam

bedtools bamtobed -i $sorted_bam > $sorted_bed

preseq c_curve -o $preseq $sorted_bed
