
#! /bin/bash
# fastq files must follow naming convention eg.
# libraryID_S1_L001_R1_001.fastq.gz
# TO RUN:
# bash cellrangerRnaCount.sh library_id path/to/fastq_dir
# eg. bash cellrangerRnaCount.sh libraryID path/to/fastq

# count the fastq files 
cellranger count \
--id=$1 \
--transcriptome="GRCh38-1.2.0_premrna" \
--fastqs=$2 \
--nosecondary


