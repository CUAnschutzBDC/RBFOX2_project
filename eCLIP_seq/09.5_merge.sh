# navigate to the output folder 08_samtools_sort

module load samtools

samtools merge ./RBFOX2_CLIP_merged.bam ./Rbfox2_CLIP_S24_L001_collapse_sorted.bam ./RBFOX2_CLIP_S31_L002_collapse_sorted.bam
samtools merge ./IgG_CLIP_merged.bam ./IgG_CLIP_S25_L001_collapse_sorted.bam ./IgG_CLIP_S30_L002_collapse_sorted.bam
samtools merge ./RBFOX2_INPUT_merged.bam ./Rbfox2_IN_S26_L001_collapse_sorted.bam ./RBFOX2_INPUT_S29_L002_collapse_sorted.bam
samtools merge ./IgG_INPUT_merged.bam ./IgG_IN_S27_L001_collapse_sorted.bam ./IgG_INPUT_S28_L002_collapse_sorted.bam

#then

samtools index ./RBFOX2_CLIP_merged.bam
samtools index ./IgG_CLIP_merged.bam
samtools index ./RBFOX2_INPUT_merged.bam
samtools index ./IgG_INPUT_merged.bam