#!/usr/bin/env bash

#BSUB -J barcode_col[1-8]
#BSUB -o logs/barcode_col_%J.out
#BSUB -e logs/barcode_col_%J.err
#BSUB -R "select[mem>64] rusage[mem=64] " 
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

python_script=collapse_barcode_pe.py

bam_file=./output/06_STAR_genome_output/${sample}.bamAligned.out.bam

out_file=./output/07_barcode_collapse/${sample}_collapse.bam

metrics_file=./output/07_barcode_collapse/${sample}_collapse.bam.metrics

echo $sample

python $python_script --bam $bam_file --out_file $out_file --metrics_file $metrics_file

# last modified on 6/30/23 for use in the replicate analysis