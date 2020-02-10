
#! /bin/bash
# this script will generate a premrna reference for counting snRNAseq libraries with hg38
# TO RUN:
# bash cellrangerRnaMkref.sh 

create custom reference that will count intronic reads
awk 'BEGIN{FS="\t"; OFS="\t"} $3 == "transcript"{ $3="exon"; print}' \
      refdata-cellranger-GRCh38-3.0.0/genes/genes.gtf > GRCh38-1.2.0.premrna.gtf

cellranger mkref \
--genome="GRCh38-1.2.0_premrna" \
--fasta="refdata-cellranger-GRCh38-3.0.0/fasta/genome.fa" \
--genes="GRCh38-1.2.0.premrna.gtf" \
--nthreads=8 \
--memgb=64



