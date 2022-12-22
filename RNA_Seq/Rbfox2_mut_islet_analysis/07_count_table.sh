#!/usr/bin/env bash

#BSUB -J featureCounts[1-8]
#BSUB -o logs/featureCounts_%J.out
#BSUB -e logs/featureCounts_%J.err
#BSUB -R "select[mem>24] rusage[mem=24] " 
#BSUB -q rna
#BSUB -n 6

module load samtools

SAMPLES=(
NT1_S5
NT3_S6
NT5_S7
NT6_S8
RbFox2_1_S1
RbFox2_3_S2
RbFox2_5_S3
RbFox2_6_S4
)

module load subread

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}


input=MIN6_RBFOX2_KD/counts/${sample}_counts.txt

echo $sample

python $python_script --file_list $input --output_file MIN6_RBFOX2_KD/MIN6_count_table.txt 

