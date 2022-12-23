#!/usr/bin/env bash

#BSUB -J featureCounts[1-10]
#BSUB -o logs/featureCounts_%J.out
#BSUB -e logs/featureCounts_%J.err
#BSUB -R "select[mem>24] rusage[mem=24] " 
#BSUB -q rna
#BSUB -n 6

module load samtools

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


module load subread

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}

GSEA=GSE183247

input=./$GSEA/${sample}_sorted.bam
output=./$GSEA/${sample}_counts.txt

echo $sample

featureCounts -a ../GEO_Datasets/ref/star/mouse/Grm38/gencode.vM25.primary_assembly.annotation.gtf --extraAttributes 'gene_name,gene_biotype' -s 2 -p -B -o $output $input