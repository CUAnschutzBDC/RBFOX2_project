#!/usr/bin/env bash

#BSUB -J rMATS
#BSUB -o logs/rMATS%J.out
#BSUB -e logs/rMATS%J.err
#BSUB -R "select[mem>16] rusage[mem=16]" 
#BSUB -q rna
#BSUB -n 6


mkdir -p logs

# Need to use python2 not python3
module load python/2.7.15

# Directory with input files
dat='/beevol/home/mossnico/Bulk_RNA-Seq'

# Path to rMATS
# rmats="$dat/rMATS-turbo-Mac-UCS2/rmats.py"
# rmats='/beevol/home/sheridanr/src/rMATS.4.0.2/rMATS-turbo-Linux-UCS2/rmats.py'
rmats='/cluster/software/modules-sw/rMATS/4.0.2/rMATS-turbo-Linux-UCS2/rmats.py'

# Lists of bam files
b1="$dat/b1_islet_KW.txt"
b2="$dat/b2_islet_KW.txt"

# GTF
gtf="../GEO_Datasets/ref/star/mouse/Grm38/gencode.vM25.primary_assembly.annotation.gtf"

# Output directory
out='Islet_Rbfox2/Islet_rMATS_output_v4.0.2'
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
    --readLength 150

