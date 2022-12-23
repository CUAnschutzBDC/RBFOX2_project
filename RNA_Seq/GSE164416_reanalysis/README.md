# GSE164416 Human Islet Bulk-RNA-Seq Analysis
This document describes the analysis pipeline for the GSE164416 dataset comparing human islets from donors with and without T2D. The scripts found in the *src* directory were used to conduct this analysis as part of study of RBFOX2 regulation.

Written by Nicole Moss

### Datasets

Human Pancreatic Islet RNA-Seq [GSE164416](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164416)

### Analysis Instructions for Initial RNA-Seq Analysis
1. fastq files obtained from GEO using [***fastq-dump***]() 
2. check quality and potentially remove adapters [***fastqc***](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and [***cutadapt***](https://cutadapt.readthedocs.io/en/stable/)
3. Align to the genome using [***STAR Aligner***](https://github.com/alexdobin/STAR)
4. Manipulate alignments using [***samtools***](http://samtools.sourceforge.net/)
5. Optional [***featureCounts***](https://subread.sourceforge.net/) and count table 
6. Alternative splicing anlysis to identify alternatively spliced genes using [***rMATS***](https://rnaseq-mats.sourceforge.net/rmats4.0.2/)

#### Downstream rMATS Analysis
Here we filter the results of the rMATS analysis for events with read counts > 20 across samples, FDR < 0.05, and |ΔPSI| > 0.01. Follow-up [kmer](https://github.com/CUAnschutzBDC/RBFOX2_project/tree/main/kmer_analysis) analysis. 


#### Downstream Differential Expression
Differential expression of genes was determined by [**DESeq2**](https://bioconductor.org/packages/release/bioc/html/DESeq2.html).



