#!/usr/bin/env bash

#BSUB -J kmer[1]
#BSUB -o logs/kmer_%J.out
#BSUB -e logs/kmer_%J.err
#BSUB -R "select[mem>64] rusage[mem=64] " 
#BSUB -q rna
#BSUB -n 6


python kmer_sequence_options_human.py -i /beevol/home/mossnico/GEO_Datasets/GSE164416_human_T2D/kmer