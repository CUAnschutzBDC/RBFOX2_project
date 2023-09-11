#!/usr/bin/env bash

#BSUB -J cut1[1-8]
#BSUB -o logs/cut_%J.out
#BSUB -e logs/cut_%J.err
#BSUB -R "select[mem>16] rusage[mem=16] " 
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
cut1_dir=output/02_adapter_cut1 
cut2_dir=output/02_adapter_cut2

r1_file=fastq/${sample}_R1_001.fastq.gz 
r2_file=fastq/${sample}_R2_001.fastq.gz 

r1_cut1=$cut1_dir/${sample}_R1_001_adapterTrim.gz 
r2_cut1=$cut1_dir/${sample}_R2_001_adapterTrim.gz 
cut1_metrics=$cut1_dir/${sample}_adapterTrim.gz.metrics 

r1_cut2=$cut2_dir/${sample}_R1_001_adapterTrim.gz 
r2_cut2=$cut2_dir/${sample}_R2_001_adapterTrim.gz 
cut2_metrics=$cut2_dir/${sample}_adapterTrim.gz.metrics 

echo $sample

mkdir $cut1_dir 
mkdir $cut2_dir

cutadapt -f fastq --match-read-wildcards --times 1 -e 0.1 -O 1 --quality-cutoff 6 -m 18 -a NNNNNAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -g CTTCCGATCTACAAGTT -g CTTCCGATCTTGGTCCT -A AACTTGTAGATCGGA -A AGGACCAAGATCGGA -A ACTTGTAGATCGGAA -A GGACCAAGATCGGAA -A CTTGTAGATCGGAAG -A GACCAAGATCGGAAG -A TTGTAGATCGGAAGA -A ACCAAGATCGGAAGA -A TGTAGATCGGAAGAG -A CCAAGATCGGAAGAG -A GTAGATCGGAAGAGC -A CAAGATCGGAAGAGC -A TAGATCGGAAGAGCG -A AAGATCGGAAGAGCG -A AGATCGGAAGAGCGT -A GATCGGAAGAGCGTC -A ATCGGAAGAGCGTCG -A TCGGAAGAGCGTCGT -A CGGAAGAGCGTCGTG -A GGAAGAGCGTCGTGT -o $r1_cut1 -p $r2_cut1 $r1_file $r2_file > $cut1_metrics

cutadapt -f fastq --match-read-wildcards --times 1 -e 0.1 -O 5 --quality-cutoff 6 -m 18 -A AACTTGTAGATCGGA -A AGGACCAAGATCGGA -A ACTTGTAGATCGGAA -A GGACCAAGATCGGAA -A CTTGTAGATCGGAAG -A GACCAAGATCGGAAG -A TTGTAGATCGGAAGA -A ACCAAGATCGGAAGA -A TGTAGATCGGAAGAG -A CCAAGATCGGAAGAG -A GTAGATCGGAAGAGC -A CAAGATCGGAAGAGC -A TAGATCGGAAGAGCG -A AAGATCGGAAGAGCG -A AGATCGGAAGAGCGT -A GATCGGAAGAGCGTC -A ATCGGAAGAGCGTCG -A TCGGAAGAGCGTCGT -A CGGAAGAGCGTCGTG -A GGAAGAGCGTCGTGT -o $r1_cut2 -p $r2_cut2 $r1_cut1 $r2_cut1 > $cut2_metrics

#last updated on 6/29/23 for use in the RBFOX2-eCLIP replicate analysis