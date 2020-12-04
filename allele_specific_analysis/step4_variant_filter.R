#!/usr/bin/env Rscript
# TO RUN INTERACTIVELY:
# export LSF_DOCKER_VOLUMES="$HOME:$HOME \
# $STORAGE1/project/vcf:$HOME/vcf \
# $STORAGE1/project/analysis/control:$HOME/outs \
# $SCRATCH1:$SCRATCH1"
# bsub -Is -G compute-parkerw -R 'rusage[mem=64GB]' -n 10 -q general-interactive -a 'docker(p4rkerw/asep:latest)' /bin/bash 

# TO RUN DETACHED:
# export LSF_DOCKER_VOLUMES="$HOME:$HOME \
# $STORAGE1/project/vcf:$HOME/vcf \
# $STORAGE1/project/analysis/control:$HOME/outs \
# $SCRATCH1:$SCRATCH1"
# SAMPLE_ARRAY=(1 2 3 4 5)
# for SAMPLE in ${SAMPLE_ARRAY[*]}
# do
# bsub -Is -G compute-parkerw -R 'rusage[mem=64GB]' -n 10 -q general-interactive -a 'docker(p4rkerw/asep:latest)' \
# Rscript $REPO/allele_specific_analysis/step4_variant.filter.R Control_${SAMPLE} atac_genotype
# done

# TO RUN LOCALLY
# docker run \
# --workdir $HOME \
# -v $HOME:$HOME \
# -v /g/project/vcf:$HOME/vcf \
# -v /g/project/analysis/control:$HOME/outs \
# -it p4rkerw/asep:latest /bin/bash

args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args) < 2) {
  stop("Please provide a sample name and output directory positional args.
    Output directory can be rna_genotype or atac_genotype.
    Example: script.R Control_1 atac_genotype", call.=FALSE)
} 

library(VariantAnnotation)
library(dplyr)
library(stringr)
library(data.table)
library(ggplot2)
library(here)

sampleName <- args[1]
output_dir <- args[2]
file <- paste0(sampleName,".variants.funcotated.vcf")
vcf <- readVcf(here("vcf", file), "hg38")

# bind all annotation into a df
df.anno <- fread(here("vcf", paste0(sampleName,".formatted.variant.csv"))) 

# remove bracket [ from gene name]
df.anno$Gencode_27_hugoSymbol <- gsub("\\[","", df.anno$Gencode_27_hugoSymbol)

# select desired anno fields
df.sel <- dplyr::select(df.anno, 
                        FILTER,
                        Gencode_27_chromosome,
                        Gencode_27_start,
                        Gencode_27_hugoSymbol, 
                        Gencode_27_variantClassification,
                        Gencode_27_codonChange,
                        gnomAD_exome_AF,
                        gnomAD_genome_AF,
                        Gencode_27_transcriptExon)

contigs <- paste0("chr", c(seq(1:22),"X","Y"))

df.sel$Gencode_27_codonChange <- na_if(df.sel$Gencode_27_codonChange, "")

# convert gnomad fields to numeric
# levels(df.sel$gnomAD_exome_AF)[1] <- NA
# df.sel$gnomAD_exome_AF <- gsub("_", " ", df.sel$gnomAD_exome_AF) %>%
#   word(1)
# df.sel$gnomAD_exome_AF[is.na(df.sel$gnomAD_exome_AF)] <- 0
# df.sel$gnomAD_exome_AF <- as.numeric(df.sel$gnomAD_exome_AF)
# df.sel$gnomAD_exome_AF[is.na(df.sel$gnomAD_exome_AF)] <- 0
# 
# levels(df.sel$gnomAD_genome_AF)[1] <- NA
# df.sel$gnomAD_genome_AF <- gsub("_", " ", df.sel$gnomAD_genome_AF) %>%
#   word(1)
# df.sel$gnomAD_genome_AF[is.na(df.sel$gnomAD_genome_AF)] <- 0
# df.sel$gnomAD_genome_AF <- as.numeric(df.sel$gnomAD_genome_AF)
# df.sel$gnomAD_exome_AF[is.na(df.sel$gnomAD_exome_AF)] <- 0
# 
# change first level of factor to NA before numeric conversion o/w it will become a 1
levels(df.sel$Gencode_27_transcriptExon)[1] <- NA
df.sel$Gencode_27_transcriptExon <- as.numeric(df.sel$Gencode_27_transcriptExon)

df.filter <- dplyr::filter(df.sel, df.sel$Gencode_27_chromosome %in% contigs)
df.filter <- dplyr::filter(df.sel, !is.na(Gencode_27_transcriptExon), FILTER == "PASS")

# hist(df.filter$gnomAD_exome_AF,
#      main="MAF of Coding Transcript Germline Variants in gnomAD Exomes",
#      xlab="Minor allele frequency (MAF)")
# 
# hist(df.filter$gnomAD_genome_AF,
#      main="MAF of Coding Transcript Germline Variants in gnomAD Genomes",
#      xlab="Minor allele frequency (MAF)")

# find the number of variants not in gnomad
nrow(dplyr::filter(df.filter, gnomAD_genome_AF == 0))

# 
# bp <- data.frame(table(df.filter$Gencode_27_variantClassification))
# ggplot(bp, aes(x=Var1,y=Freq,fill=Var1)) + geom_bar(stat="identity") +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# filter for exonic variants
vcf.mrna <- vcf[which(!is.na(df.sel$Gencode_27_transcriptExon))]

# filter for variants that PASS
vcf.mrna <- vcf.mrna[which(rowRanges(vcf.mrna)$FILTER == "PASS")]

# write the vcf to file
writeVcf(vcf.mrna, here("vcf", paste0(sampleName,".mrna.vcf")))

# filter for variants that overlap with peaks
peaks.gr <- fread(here("outs","atac_aggr_prep","step6_peaks.bed")) %>%
  dplyr::rename("chr" = V1, "start" = V2, "end" = V3) %>%
  makeGRangesFromDataFrame()

overlap <- findOverlaps(rowRanges(vcf), peaks.gr)
vcf.peaks <- vcf[queryHits(overlap)]
vcf.peaks <- vcf.peaks[which(rowRanges(vcf.peaks)$FILTER == "PASS")]

# write the vcf to file
writeVcf(vcf.peaks, here("vcf", paste0(sampleName,".peaks.vcf")))

# filter for intronic and exonic variants
vcf.premrna <- vcf[which(!is.na(df.sel$Gencode_27_transcriptExon) | df.sel$Gencode_27_variantClassification == "INTRON"),]
vcf.premrna <- vcf.premrna[which(rowRanges(vcf.premrna)$FILTER == "PASS")]
writeVcf(vcf.premrna, here("vcf", output_dir, paste0(sampleName,".premrna.vcf")))


