# eCLIP-Seq Analysis Pipeline
Adapted from Van Nostrand, et al., 2016

### General Analysis Pipeline

* [**01_fastqc1.sh**](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) (v0.11.9)
* **02_cutadapt.sh**
* **03_STAR_rmRep.sh**
* [**04_fastqc2.sh**](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) (v0.11.9)
* **05_fastq_sort.sh**
* [**06_STAR_genome.sh**](https://github.com/alexdobin/STAR) (2.7.9a)
* [*06b_preseq.sh*](http://smithlabresearch.org/software/preseq/) - Optional assessment of sequencing satruation. 
* **07_collapse_barcode.sh** - Using *collapse_barcode_pe.py* reqritten in python. 
* [**08_samtools_sort.sh**](http://www.htslib.org/) (Li et al., 2009) 
* [**09_samtools_index.sh**](http://www.htslib.org/) (Li et al., 2009)
* **09.5_merge.sh** Optional combine replicate sequencing experiments. 
* **10_samtools_view.sh** OR **10_samtools_view_merge.sh**
* **11_samtools_index.sh** OR **11_samtools_index_merge.sh**
* **12_bigwig2.sh** OR **12_bigwig_merge.sh**
* **13_clipper2_.sh** OR **13_clipper2_merge.sh** (https://github.com/YeoLab/clipper)
* **14_SMI.sh** Peak normalization has been modified from the original publication to call peaks based on the appropriate read for the Eclipse BioInnovations RBP-eCLIP and library prep kit [(https://eclipsebio.com/services-kits/rbp-eclip/)]. Details on this code can be found at https://github.com/kwells4/eclip_normalization OR DOI: 10.5281/zenodo.8335301.

Replicate comparison by IDR (https://github.com/kwells4/eclip_idr OR DOI: 10.5281/zenodo.8335366)
