# Bulk RNA-Seq Analysis Pipeline

### Experiments
**Rbfox2-KD in MIN6 Cells**
  Experimental Design
  Rbfox2 KD by pooled siRNA in MIN6 compared to pool of 4 non-target control siRNAs. 
\
\
####Rbfox2-mut in mouse islets**
Bulk RNA-Seq comparing islets isolated from Pdx1:Cre;Rbfox2fl/lf (experimental) 8-week mice with Rbfox2fl/fl (control) 8-week mice.  
|   *822* - Male Rbfox2fl/fl 
|   *824* - Female Rbfox2fl/fl
|   *825* - Female Pdx1:Cre;Rbfox2fl/fl
|   *826* - Male Rbfox2fl/fl
|   *827* - Male Pdx1:Cre;Rbfox2fl/fl
|   *828* - Male Pdx1:Cre;Rbfox2fl/fl
\
\
### Sequencing Parameters
| CU Anschutz Genomics and Micorarray Core
| Poly A Selected Total RNA - NovaSEQ 6000 (Paired End 150 cycle 2x150)
\
\
### Analysis Pipeline
  FastQC (v0.11.9) 
  Cutadapt (Martin, 2011)
  STAR aligner (2.7.9a)
  samtools (Li et al., 2009)
  featureCounts (subread 1.6.2) (Liao et al., 2014)
  DESeq2 (v1.34.0) 
  rMATS (v4.0.2) 
  \
**rMATS v4.0.2 Splicing Analysis**
  \
|   **b1_MIN6.txt** NT1_S5_sorted.bam,NT3_S6_sorted.bam,NT5_S7_sorted.bam,NT6_S8_sorted.bam
|   **b2_MIN6.txt** RbFox2_1_S1_sorted.bam,RbFox2_3_S2_sorted.bam,RbFox2_5_S3_sorted.bam,RbFox2_6_S4_sorted.bam
\
|   **b1.txt** 822_S34_sorted.bam,824_S35_sorted.bam,826_S37_sorted.bam
|   **b2.txt** 825_S36_sorted.bam,827_S38_sorted.bam,828_S39_sorted.bam
\
  **rMATS** generates several output files, for this analysis we will use the 'JC.txt' files that use only splice junction spanning reads to quantify alternative splicing events. The files are listed below...
\
|   **A3SS.MATS.JC.txt** - Alternative 3' Splice Site
|   **A5SS.MATS.JC.txt** - Alternative 5' Splice Site
|   **RI.MATS.JC.txt** - Retained Intron
|   **MXE.MATS.JC.txt** - Mutually Exclusive Exon
|   **SE.MATS.JC.txt** - Skipped Exon
\
  **CITATIONS**:
Shen S., Park JW., Lu ZX., Lin L., Henry MD., Wu YN., Zhou Q., Xing Y. rMATS: Robust and Flexible Detection of Differential Alternative Splicing from Replicate RNA-Seq Data. PNAS, 111(51):E5593-601. doi: 10.1073/pnas.1419161111 

Park JW., Tokheim C., Shen S., Xing Y. Identifying differential alternative splicing events from RNA sequencing data using RNASeq-MATS. Methods in Molecular Biology: Deep Sequencing Data Analysis, 2013;1038:171-179 doi: 10.1007/978-1-62703-514-9_10 

Shen S., Park JW., Huang J., Dittmar KA., Lu ZX., Zhou Q., Carstens RP., Xing Y. MATS: A Bayesian Framework for Flexible Detection of Differential Alternative Splicing from RNA-Seq Data. Nucleic Acids Research, 2012;40(8):e61 doi: 10.1093/nar/gkr1291
