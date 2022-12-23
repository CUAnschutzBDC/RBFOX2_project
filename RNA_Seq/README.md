# Bulk RNA-Seq Analysis Pipeline

### General Analysis Pipeline
* [*fastq-dump*](https://rnnh.github.io/bioinfo-notebook/docs/fastq-dump.html)
* [*FastQC*](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) (v0.11.9) 
* [*cutadapt*](https://cutadapt.readthedocs.io/en/stable/) (Martin, 2011)
* [STAR aligner](https://github.com/alexdobin/STAR) (2.7.9a)
* [*samtools*](http://www.htslib.org/) (Li et al., 2009)
* [*featureCounts*](https://academic.oup.com/bioinformatics/article/30/7/923/232889) (subread 1.6.2) (Liao et al., 2014)
* [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) (v1.34.0) 
* [rMATS](https://rnaseq-mats.sourceforge.net/rmats4.0.2/) (v4.0.2) 
 
### Experiments

[**GSE164416_reanalysis**](https://github.com/CUAnschutzBDC/RBFOX2_project/tree/main/RNA_Seq/GSE164416_reanalysis) - This experiment reanalyzes data from [GSE164416](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164416) which contains bulk islet RNA-Seq from heathy and diabetic donors. This dataset was used to identify *RBFOX2* expression changes in T2D and the enrichment of RBFOX2 binding sequence (GCAUG) near alternative exons. 

[**GSE183247_reanalysis**](https://github.com/CUAnschutzBDC/RBFOX2_project/tree/main/RNA_Seq/GSE183247_reanalysis) - This experiment reanalyzes data from [GSE183247](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE183247) which contains bulk islet RNA-Seq from obese diabetic mice (NZO) and obese non-diabetic mice (*ob/ob*). This dataset was used to identify *Rbfox2* expression changes in T2D and the enrichment of RBFOX2 binding sequence (GCAUG) near alternative exons. 

[**HPAP_sorted_beta_reanalysis**](https://github.com/CUAnschutzBDC/RBFOX2_project/tree/main/RNA_Seq/HPAP_sorted_beta_reanalysis) - This experiment reanalyzes data from [PANC-DB](https://hpap.pmacs.upenn.edu/). Here we use the sorted human β cells to identify changes in *RBFOX2*.

[**Rbfox2_KD_MIN6_analysis**](https://github.com/CUAnschutzBDC/RBFOX2_project/tree/main/RNA_Seq/Rbfox2_KD_MIN6_analysis) - This experiment is specific to our study of the role of RBFOX2 in pancreatic β cells and will be available at GSE221277.

* In this analysis we compare mouse insulinoma (MIN6) cells treated with a pool of non-target siRNA (NT) or *Rbfox2* directed siRNA (KD) for 48h. 

* PolyA Selected Total RNA, NovaSeq 6000 150bp Paired End, ~40 million Reads

[**Rbfox2_mut_islet_analysis**](https://github.com/CUAnschutzBDC/RBFOX2_project/tree/main/RNA_Seq/Rbfox2_mut_islet_analysis) - This experiment is specific to our study of the role of RBFOX2 in mouse pancreas and will be available at GSE221277.

* In this analysis we compare bulk RNA-Seq from *Rbfox2*-mut islets (Pdx1:CRE; Rbfox2 fl/fl) to the corresponding non-CRE controls (Rbfox2 fl/fl). 

* PolyA Selected Total RNA, NovaSeq 6000 150bp Paired End, ~40 million Reads
