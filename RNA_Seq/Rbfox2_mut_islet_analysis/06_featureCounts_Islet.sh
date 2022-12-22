#!/usr/bin/env bash

#BSUB -J featureCounts[1-6]
#BSUB -o logs/featureCounts_%J.out
#BSUB -e logs/featureCounts_%J.err
#BSUB -R "select[mem>24] rusage[mem=24] " 
#BSUB -q rna
#BSUB -n 6

module load samtools

module load samtools

SAMPLES=(
822_S34
824_S35
825_S36
826_S37
827_S38
828_S39
)

module load subread

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}


input=${sample}_sorted.bam
output=${sample}_counts.txt

echo $sample

featureCounts -a ../GEO_Datasets/ref/star/mouse/Grm38/gencode.vM25.primary_assembly.annotation.gtf --extraAttributes 'gene_name,gene_biotype' -s 2 -p -B -o $output $input