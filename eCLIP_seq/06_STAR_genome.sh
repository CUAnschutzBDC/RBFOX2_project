#!/usr/bin/env bash

#BSUB -J STAR_genome[1-8]
#BSUB -o logs/STAR_genome_%J.out
#BSUB -e logs/STAR_genome_%J.err
#BSUB -R "select[mem>64] rusage[mem=64] " 
#BSUB -q rna
#BSUB -n 6

module load STAR/2.7.9a

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
STAR_dir=./output/06_STAR_genome_output

r1_file=./output/05_unmapped_sorted_fastq/${sample}_rep.bamUnmapped.out.mate1_sorted
r2_file=./output/05_unmapped_sorted_fastq/${sample}_rep.bamUnmapped.out.mate2_sorted

STAR_file=$STAR_dir/${sample}.bam

echo $sample

STAR --runMode alignReads --runThreadN 16 --genomeDir ../ref/star/mouse/Grm38/index --genomeLoad LoadAndRemove --readFilesIn $r1_file $r2_file --outSAMunmapped Within --outFilterMultimapNmax 1 --outFilterMultimapScoreRange 1 --outFileNamePrefix $STAR_file --outSAMattributes All --outSAMtype BAM Unsorted --outFilterType BySJout --outReadsUnmapped Fastx --outFilterScoreMin 10 --outSAMattrRGline ID:foo --alignEndsType EndToEnd 

# last modified on 6/30/23 for use in the replicate analysis 
