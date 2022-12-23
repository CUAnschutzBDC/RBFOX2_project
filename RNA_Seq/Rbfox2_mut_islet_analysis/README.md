# Rbfox2-mut Mouse Islets
The goal of this experiment is to establish a model to investigate RBFOX2 sensitive alternative splicing in mouse islets. 

Written by Nicole Moss

### Datasets

This experiment is specific to our study of the role of RBFOX2 in pancreatic β cells and will be available at GSE221277.

* In this analysis we compare bulk-RNA-Seq from 8-week Rbfox2_mut islets (Pdx1:CRE;Rbfox2 fl/fl) to the coresponding non-cre control (Rbfox2 fl/fl), n = 3 for each condition.

* PolyA Selected Total RNA, NovaSeq 6000 150bp Paired End, ~40 million Reads


### Analysis Instructions for Initial RNA-Seq Analysis
1. check quality [***fastqc***](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
2. remove adapters [***cutadapt***](https://cutadapt.readthedocs.io/en/stable/)
3. re-check quality [***fastqc***](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
4. Align to the genome using [***STAR Aligner***](https://github.com/alexdobin/STAR)
5. Manipulate alignments using [***samtools***](http://samtools.sourceforge.net/)
6. Optional [***featureCounts***](https://subread.sourceforge.net/) and count table 7. generate count table
8. Alternative splicing anlysis to identify alternatively spliced genes using [***rMATS***](https://rnaseq-mats.sourceforge.net/rmats4.0.2/)

#### Downstream rMATS Analysis
Here we filter the results of the rMATS analysis for events with read counts > 20 across samples, FDR < 0.05, and |ΔPSI| > 0.01. Follow-up [kmer](https://github.com/CUAnschutzBDC/RBFOX2_project/tree/main/kmer_analysis) analysis. 


#### Downstream Differential Expression
Differential expression of genes was determined by [**DESeq2**](https://bioconductor.org/packages/release/bioc/html/DESeq2.html).


