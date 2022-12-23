# eCLIP Peak Sliding Window Analysis
This function will count overlapping peaks within a set distance from exon boundaries. 

Written by Nicole Moss and Kristen Wells

### Datasets 
Initial datasets for this analysis are the normalized and IgG controlled output of the eCLIP analysis pipeline and the filterd output of rMATS casette exon (or skipped exon) and mutually exclusive exon. 

### Function

**peak_exon_overlap**(rMATS_data, eCLIP_data, nt_start, nt_end, by)

* [*rMATS_data*](https://github.com/CUAnschutzBDC/RBFOX2_project/tree/main/RNA_Seq/Rbfox2_KD_MIN6_analysis)
* [*eCLIP_data*](https://github.com/CUAnschutzBDC/RBFOX2_project/tree/main/eCLIP_seq)
* *nt_start* to *nt_end* : the range that you would like to investigate the overlap of peaks
* *by* : window size within the range that you would like to consider
  
  
