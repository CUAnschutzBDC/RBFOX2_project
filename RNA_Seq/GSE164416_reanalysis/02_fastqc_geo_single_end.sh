#!/usr/bin/env bash

#BSUB -J fastqc1[1-73]
#BSUB -o logs/fastqc2_%J.out
#BSUB -e logs/fastqc2_%J.err
#BSUB -R "select[mem>40] rusage[mem=40] " 
#BSUB -q rna
#BSUB -n 6

module load fastqc/0.11.7

SAMPLES=(
SRR12885579
SRR12885580
SRR12885581
SRR12885582
SRR12885583
SRR12885584
SRR12885585
SRR12885586
SRR12885587
SRR12885588
SRR12885589
SRR12885590
SRR12885591
SRR12885592
SRR12885593
SRR12885594
SRR12885595
SRR12885596
SRR12885597
SRR12885598
SRR12885599
SRR12885600
SRR12885601
SRR12885602
SRR12885603
SRR12885604
SRR12885605
SRR12885606
SRR12885607
SRR12885608
SRR12885609
SRR12885610
SRR12885611
SRR12885612
SRR12885613
SRR12885614
SRR12885615
SRR12885616
SRR12885617
SRR12885618
SRR12885619
SRR12885620
SRR12885621
SRR12885622
SRR12885623
SRR12885624
SRR12885625
SRR12885626
SRR12885627
SRR12885628
SRR12885629
SRR12885630
SRR12885631
SRR12885632
SRR12885633
SRR12885634
SRR12885635
SRR12885636
SRR12885637
SRR12885638
SRR12885639
SRR12885640
SRR12885641
SRR12885642
SRR12885643
SRR12885644
SRR12885645
SRR12885646
SRR12885647
SRR12885648
SRR12885649
SRR12885650
SRR12885651
)


sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}

echo $sample

GSEA=GSE159984

r1_file=$GSEA/${sample}_1.fastq

fastqc $r1_file -o fastqc

