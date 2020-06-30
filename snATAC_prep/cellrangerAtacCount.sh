#! /bin/bash
# fastq files must follow naming convention eg.
# libraryID_S1_L001_R1_001.fastq.gz
# TO RUN:
# bash cellrangerAtacCount.sh library_id path/to/fastq_dir
# eg. bash cellrangerAtacCount.sh libraryID path/to/fastq

# count the fastq files 
# cellranger-atac count \
# --id=ADN2_034 \
# --reference=/home/parkerw/reference/refdata-cellranger-atac-GRCh38-1.2.0 \
# --fastqs=. \
# --sample=ADN2_034

# count the fastq files 
cellranger-atac count \
--id=$1 \
--reference=/home/parkerw/reference/refdata-cellranger-atac-GRCh38-1.2.0 \
--fastqs=$2 \
--sample=$1

