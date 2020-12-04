#! /bin/bash
# fastq files must follow naming convention eg.
# libraryID_S1_L001_R1_001.fastq.gz

# count the fastq files 
cellranger-atac count \
--id=$1 \
--reference=$HOME/reference/refdata-cellranger-atac-GRCh38-1.2.0 \
--fastqs=$2 \
--sample=$1

