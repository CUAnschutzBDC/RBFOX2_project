# kmer Analysis
This document describes the kmer analysis used to identify kmer enrichment near alternative exons. The scripts found in the *src* directory were used to conduct this analysis as part of study of RBFOX2 regulation in both mouse and human type 2 diabetes (T2D). 

Written by Nicole Moss

### Datasets

Human Pancreatic Islet RNA-Seq [GSE164416](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164416)

Mouse Pancreatic Islet RNA-Seq [GSE183247](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE183247)

### Analysis Instructions for Initial RNA-Seq Analysis 
(from the **RNA_Seq** repo)\

1. fastq files obtained from GEO using ***fastq-dump*** 
2. check quality and potentially remove adapters ***fastqc*** and ***cutadapt***
3. Align to the genome using ***STAR Aligner***
4-5. Manipulate alignments using ***samtools***
6-7. Optional ***featureCounts*** and count table 
8. Alternative splicing anlysis to identify alternatively spliced genes using ***rMATS***

### kmer Analysis Sepcific Processing 
(referenced srcipts are in the *src* directory)

1. Starting with the skipped exon output files from rMATS (*SE.MATS.JCEC.txt* or *SE.MATS.JC.txt*), generate a list of coordinates for alternative exons that are skipped (FDR < 0.05 and IncLevelDifference < -0.01), included (FDR < 0.05 and IncLevelDifference > 0.01), or insensitive (FDR > 0.05). ***kmer_analysis_coordinates_human.R*** or ***kmer_analysis_coordinates_mouse.R*** this will generate .txt files with coordinates of each exon.
2. Translate coordinates into sequences and tally the occurances of every possible kmer (5-mer in this case) within a given distance from the exon boundary (here 200nt) using the ***kmer_sequence_options_human.py*** or ***kmer_sequence_options_mouse.py*** returning a fastq sequence file and kmer_coults.txt file for each condition.
3. Determine enrichment of kmer counts relative to insensitive exons using ***kmer_enrichment_plot_mouse.R*** or ***kmer_enrichment_plot_mouse.R***

