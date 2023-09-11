#!/usr/bin/env bash

#BSUB -J SMI[1]
#BSUB -o logs/SMI2_%J.out
#BSUB -e logs/SMI2_%J.err
#BSUB -R "select[mem>64] rusage[mem=64] " 
#BSUB -q rna
#BSUB -n 6

module load samtools
module load R

python peak_input_normalization.py -m manfest_file.txt -o output/14_SMI_results/merged 
