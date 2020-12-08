
#! /bin/bash
# fastq files must follow naming convention eg.
# libraryID_S1_L001_R1_001.fastq.gz

# count the fastq files 
cellranger count \
--id=$1 \
--transcriptome="GRCh38-1.2.0_premrna" \
--fastqs=$2 \
--nosecondary


