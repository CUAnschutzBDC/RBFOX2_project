#!/usr/bin/env bash

#BSUB -J featureCounts[1-3]
#BSUB -o logs/featureCounts_%J.out
#BSUB -e logs/featureCounts_%J.err
#BSUB -R "select[mem>24] rusage[mem=24] " 
#BSUB -q rna
#BSUB -n 6

module load samtools

SAMPLES=(
HPAP-051
HPAP-053
HPAP-057
)

module load subread

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}

GSEA=HPAP

input=./$GSEA/${sample}_sorted.bam
output=./$GSEA/${sample}_counts.txt

echo $sample

featureCounts -a ./ref/star/human/GRCh38/gencode.v41.annotation.gtf --extraAttributes 'gene_name,gene_biotype' -s 2 -p -B -o $output $input