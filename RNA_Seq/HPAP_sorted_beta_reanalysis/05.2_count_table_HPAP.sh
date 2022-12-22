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


input=HPAP/${sample}_counts.txt

echo $sample

python $python_script --file_list $input --output_file HPAP/count_table.txt 