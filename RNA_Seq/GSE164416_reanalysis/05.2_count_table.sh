#!/usr/bin/env bash

#BSUB -J featureCounts[1-57]
#BSUB -o logs/featureCounts_%J.out
#BSUB -e logs/featureCounts_%J.err
#BSUB -R "select[mem>24] rusage[mem=24] " 
#BSUB -q rna
#BSUB -n 6

module load samtools

SAMPLES=(
SRR13380421
SRR13380425
SRR13380426
SRR13380427
SRR13380430
SRR13380431
SRR13380434
SRR13380435
SRR13380438
SRR13380441
SRR13380443
SRR13380445
SRR13380446
SRR13380448
SRR13380452
SRR13380457
SRR13380458
SRR13380466
SRR13380470
SRR13380471
SRR13380472
SRR13380477
SRR13380479
SRR13380480
SRR13380481
SRR13380482
SRR13380483
SRR13380484
SRR13380486
SRR13380488
SRR13380490
SRR13380495
SRR13380497
SRR13380499
SRR13380500
SRR13380501
SRR13380502
SRR13380503
SRR13380506
SRR13380508
SRR13380511
SRR13380512
SRR13380513
SRR13380516
SRR13380518
SRR13380520
SRR13380521
SRR13380525
SRR13380528
SRR13380529
SRR13380531
SRR13380533
SRR13380535
SRR13380537
SRR13380538
SRR13380539
SRR13380541
)7
)

module load subread

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}


input=${sample}_counts.txt

echo $sample

python $python_script --file_list $input --output_file /count_table.txt 