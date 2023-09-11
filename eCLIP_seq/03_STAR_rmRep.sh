#!/usr/bin/env bash

#BSUB -J STAR_rmRep[1-8]
#BSUB -o logs/STAR_rmRep_%J.out
#BSUB -e logs/STAR_rmRep_%J.err
#BSUB -R "select[mem>16] rusage[mem=16] " 
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
rmRep_dir=output/03_rmRep

r1_file=./output/02_adapter_cut2/${sample}_R1_001_adapterTrim.gz
r2_file=./output/02_adapter_cut2/${sample}_R2_001_adapterTrim.gz

rmRep_file=$rmRep_dir/${sample}_rep.bam

mkdir $rmRep_dir
echo $sample


STAR --runMode alignReads --runThreadN 16 --genomeDir ../ref/star/mouse/Grm38/repeats/STAR_rep_index --genomeLoad LoadAndRemove --readFilesIn $r1_file $r2_file --outSAMunmapped Within --outFilterMultimapNmax 30 --outFilterMultimapScoreRange 1 --outFileNamePrefix $rmRep_file --outSAMattributes All --readFilesCommand zcat --outStd BAM_Unsorted --outSAMtype BAM Unsorted --outFilterType BySJout --outReadsUnmapped Fastx --outFilterScoreMin 10 --outSAMattrRGline ID:foo --alignEndsType EndToEnd 


#last updated on 6/29/23 for use in the RBFOX2-eCLIP replicate analysis