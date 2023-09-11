#!/usr/bin/env bash

#BSUB -J make_bigwig[1-4]
#BSUB -o logs/bigwig_%J.out
#BSUB -e logs/bigwig_%J.err
#BSUB -R "select[mem>64] rusage[mem=64] " 
#BSUB -q rna
#BSUB -n 2

SAMPLES=(
IgG_CLIP
IgG_INPUT
RBFOX2_CLIP
RBFOX2_INPUT
)

module load samtools

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}

r1_file=./output/08_samtools_sort/${sample}_merged.bam

out=./output/12_bigwig/${sample}_merged.bw

echo $sample

samtools index $r1_file

bamCoverage -b $r1_file \
		-o $out \
		-p max/2 \
		--binSize 6 \
		--effectiveGenomeSize 2652783500 \
		--normalizeUsing RPKM

# last modified on 7/7/23 for use in merged replicate analysis
