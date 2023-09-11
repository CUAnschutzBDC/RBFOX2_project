#!/usr/bin/env bash

#BSUB -J make_bigwig[1-8]
#BSUB -o logs/bigwig_%J.out
#BSUB -e logs/bigwig_%J.err
#BSUB -R "select[mem>64] rusage[mem=64] " 
#BSUB -q rna
#BSUB -n 2

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
module load samtools

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}

r1_file=./output/08_samtools_sort/${sample}_collapse_sorted.bam

out=./output/12_bigwig/${sample}.bw

echo $sample

samtools index $r1_file

bamCoverage -b $r1_file \
		-o $out \
		-p max/2 \
		--binSize 6 \
		--effectiveGenomeSize 2652783500 \
		--normalizeUsing RPKM

# last modified on 7/5/23 for use in replicate analysis
