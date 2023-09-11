#!/usr/bin/env bash

#BSUB -J fastq_sort[1-8]
#BSUB -o logs/fastq_sort_%J.out
#BSUB -e logs/fastq_sort_%J.err
#BSUB -R "select[mem>40] rusage[mem=40] " 
#BSUB -q rna
#BSUB -n 6

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

r1_file=./output/03_rmRep/${sample}_rep.bamUnmapped.out.mate1
r2_file=./output/03_rmRep/${sample}_rep.bamUnmapped.out.mate2

sorted_file1=./output/05_unmapped_sorted_fastq/${sample}_rep.bamUnmapped.out.mate1_sorted
sorted_file2=./output/05_unmapped_sorted_fastq/${sample}_rep.bamUnmapped.out.mate2_sorted

echo $sample

fastq-sort --id $r1_file > $sorted_file1

fastq-sort --id $r2_file > $sorted_file2

# last modified on 6/29/23 for use in replicate analysis
# use the "alignment_counting" conda enviroment 
