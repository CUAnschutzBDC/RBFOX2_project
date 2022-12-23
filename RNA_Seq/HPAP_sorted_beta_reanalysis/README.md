# HPAP Human sorted β cells

Written by Nicole Moss

### Datasets

[**HPAP_sorted_beta_reanalysis**](https://github.com/CUAnschutzBDC/RBFOX2_project/tree/main/RNA_Seq/HPAP_sorted_beta_reanalysis) - This experiment reanalyzes data from [PANC-DB](https://hpap.pmacs.upenn.edu/). Here we use the sorted human β cells to identify changes in *RBFOX2*.

### Analysis Instructions for Initial RNA-Seq Analysis
1. fastq files obtained from GEO using [***fastq-dump***]() 
2. check quality and potentially remove adapters [***fastqc***](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and [***cutadapt***](https://cutadapt.readthedocs.io/en/stable/)
3. Align to the genome using [***STAR Aligner***](https://github.com/alexdobin/STAR)
4. Manipulate alignments using [***samtools***](http://samtools.sourceforge.net/)
5. Optional [***featureCounts***](https://subread.sourceforge.net/) and count table 

#### Limited downstream analysis due to low numbner of replicates. 

