#!/usr/bin/env bash

#BSUB -J rMATS
#BSUB -o logs/rMATS%J.out
#BSUB -e logs/rMATS%J.err
#BSUB -R "select[mem>64] rusage[mem=64]" 
#BSUB -q rna
#BSUB -n 6


mkdir -p logs

# Need to use python2 not python3
module load python/2.7.15

# Directory with input files
dat='/beevol/home/mossnico/GEO_Datasets'

# Path to rMATS
# rmats="$dat/rMATS-turbo-Mac-UCS2/rmats.py"
# rmats='/beevol/home/sheridanr/src/rMATS.4.0.2/rMATS-turbo-Linux-UCS2/rmats.py'
rmats='/cluster/software/modules-sw/rMATS/4.0.2/rMATS-turbo-Linux-UCS2/rmats.py'

# Lists of bam files
b1="$dat/b1_GSE183247.txt"
b2="$dat/b2_GSE183247.txt"

# GTF
gtf="$dat/ref/star/mouse/Grm38/gencode.vM25.primary_assembly.annotation.gtf"

# Output directory
out='rMATS_v4.0.2_output_GSE183247'
out=$(realpath "$out")

# Run rMATS
# specify number of threads with --nthread
cd "$dat"

python "$rmats" \
    --b1 "$b1" \
    --b2 "$b2" \
    --gtf "$gtf" \
    --od "$out" \
    -t paired \
    --nthread 6 \
    --readLength 125
