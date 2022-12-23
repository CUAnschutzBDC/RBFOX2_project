#!/usr/bin/env bash

#BSUB -J GEO[1-10]
#BSUB -o logs/GEO_%J.out
#BSUB -e logs/GEO_%J.err
#BSUB -R "select[mem>64] rusage[mem=64] " 
#BSUB -q rna
#BSUB -n 6

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

module load sratoolkit


sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}

echo $sample

fastq-dump --split-files $sample