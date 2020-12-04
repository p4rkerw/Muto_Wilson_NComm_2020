#!/usr/bin/env Rscript
# TO RUN INTERACTIVELY:
# export LSF_DOCKER_VOLUMES="$HOME:$HOME \
# $STORAGE1/project/analysis/control/wasp_snRNA_joint_premrna:$HOME/allele_spec \
# $SCRATCH1:$SCRATCH1"
# bsub -Is -G compute-parkerw -R 'rusage[mem=128GB]' -q general-interactive -a 'docker(p4rkerw/asep:latest)' /bin/bash 

# # TO RUN DETACHED:
# export LSF_DOCKER_VOLUMES="$HOME:$HOME \
# $STORAGE1/project/analysis/control/wasp_snRNA_joint_premrna:$HOME/allele_spec \
# $SCRATCH1:$SCRATCH1"
# bsub -G compute-parkerw -R 'rusage[mem=128GB]' -J "asep_premrna_DCTPC" -q general -o asep_premrna_DCTPC.out -a 'docker(p4rkerw/asep:latest)' \
# Rscript healthy_dev/gatk/rna_allele_specific/step7_asep_analyze.R DCTPC premrna 3
# bsub -G compute-parkerw -R 'rusage[mem=128GB]' -J "asep_mrna_DCTPC" -q general -o asep_mrna_DCTPC.out -a 'docker(p4rkerw/asep:latest)' \
# Rscript healthy_dev/gatk/rna_allele_specific/step7_asep_analyze.R DCTPC mrna 3


# RUN LOCAL INTERACTIVE:
# docker run \
# --workdir $HOME \
# -v $HOME:$HOME \
# -v /g/project/vcf/atac_genotype:$HOME/vcf \
# -v /g/project/analysis/control/wasp_snRNA_joint_premrna:$HOME/allele_spec \
# -v /g/project/analysis/control/barcodes:$HOME/barcodes \
# -v /g/project/cellranger_atac_counts/version_1.2:$HOME/counts \
# -v /g/reference/refdata-cellranger-atac-GRch38-1.2.0:$HOME/ref \
# -v /g/scratch:/scratch \
# -v $SCRATCH1:$SCRATCH1 \
# -e SCRATCH1="/g/scratch" \
# --rm -it p4rkerw/asep:latest

# RUN LOCAL INTERACTIVE RSTUDIO VIA CHROME localhost:8787
# docker run \
# -v $HOME:$HOME \
# -v /g/project:/home/rstudio/ \
# -v /g/project/vcf:/home/rstudio/vcf \
# -v /g/project/analysis/control/wasp_snRNA_joint_premrna:/home/rstudio/allele_spec \
# -v /g/project/analysis/control/barcodes:/home/rstudio/barcodes \
# -v /g/project/cellranger_atac_counts/version_1.2:/home/rstudio/counts \
# -v /g/reference/refdata-cellranger-atac-GRch38-1.2.0:/home/rstudio/ref \
# -v /g/scratch:/home/rstudio/scratch \
# --rm -p 8787:8787 -e PASSWORD=password p4rkerw/asep:latest

library(ASEP)
library(data.table)
library(dplyr)
library(stringr)
library(GenomicRanges)
library(EnsDb.Hsapiens.v86)
library(VariantAnnotation)
library(here)
library(plyranges)
library(tibble)
library(multtest)

args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args) < 2) {
  stop("Please provide celltype, vcf region, and resampling exponent as positional arguments
    Example: script.R ALLCELLS mrna 4
    Example: script.R PT premrna 3", call.=FALSE)
} 

celltype <- args[1]
vcf_region <- args[2]
resample <- args[3]

# read in the allele-specific expression read counts from snRNA ASEReadCounter following GATK identification
# of biallelic SNVs from snATAC data
celltype_ase_table <- paste0("ase.",vcf_region,".",celltype,".output.table")
files <- list.files(here("allele_spec"), pattern=celltype_ase_table, recursive=TRUE, full.names=TRUE)
files <- files[grepl("Control", files)] 
ase.ls <- lapply(files, function(x) {fread(x)})

# load the ensembl annotations
gene.coords <- genes(EnsDb.Hsapiens.v86, filter = ~ gene_biotype == "protein_coding", columns = "gene_name")
seqlevelsStyle(gene.coords) <- 'UCSC'
gene.coords <- keepStandardChromosomes(gene.coords, pruning.mode = 'coarse')

# load the snv from the count table into a GRanges object
snv.df <- lapply(seq(ase.ls), function(x){
  ase <- ase.ls[[x]]
  snv <- GRanges(
  seqnames=ase$contig,
  ranges=IRanges(start=ase$position, width=1),
  allele1=ase$refAllele,
  allele2=ase$altAllele,
  allele1_counts=ase$refCount,
  allele2_counts=ase$altCount,
  totalCount=ase$totalCount,
  id=x,
  group="Control"
  )
  # find the genes that correspond to each snv (some snv overlap with multiple genes)
  snv <- join_overlap_intersect(snv, gene.coords)
  snv$aseID <- snv$gene_name
  return(as.data.frame(snv))
}) %>% bind_rows()

# filter out positions where the minor allele count is < 5 or
# the total read count is less than 20 or
# the minor allele count is less than 5% of total read count
snv.df$minor_count <- with(snv.df, pmin(allele1_counts, allele2_counts))
snv.df <- dplyr::mutate(snv.df, prop_minor_count = minor_count/totalCount)
snv.df <- dplyr::filter(snv.df, totalCount > 20,
                        allele1_counts > 5, 
                        allele2_counts > 5,
                        prop_minor_count > 0.05)


asep.df <- dplyr::mutate(snv.df, snp=paste0(seqnames,":",start)) %>%
  dplyr::rename(gene=gene_name) %>%
  dplyr::mutate(ref=allele1_counts) %>%
  dplyr::mutate(total=allele1_counts + allele2_counts) %>%
  dplyr::select(gene, id, group, snp, ref, total) %>%
  dplyr::mutate(ref_condition=group) %>%
  dplyr::arrange(gene)

# filter the df for genes that are represented by at least 3 individuals 
ase.genes <- unique(asep.df$gene)
lookup <- sapply(ase.genes, function(lookup) {
  tmp <- asep.df[asep.df$gene == lookup,]
  num=length(unique(tmp$id))
  }, simplify=TRUE) %>%
as.data.frame() %>%
rownames_to_column(var="gene") %>%
dplyr::filter(. > 2) 
genes <- lookup$gene

# grab genes to test
asep.sub <- asep.df[asep.df$gene %in% genes, ]

# identify genes with allele specific expression across individuals
# and genes where at least 3 individuals have a snv
Sys.time()
res.df <- lapply(genes, function(gene) {
  print(gene)
  df <- asep.sub[asep.sub$gene == gene,]
  res <- tryCatch(ASE_detection(df, phased=FALSE, adaptive=TRUE,
                         n_resample=10^as.numeric(resample), parallel=TRUE, save_out=FALSE),
           error = function(e) NULL)
  closeAllConnections()
  if(!is.null(res)) {
    df.return <- data.frame(gene=res[1], pval=res[2])
    print(df.return)
    return(df.return)
  } else {
    return(NULL)
  }
}) %>% bind_rows()
Sys.time()
dir.create(here("allele_spec","results"))

# perform multiple testing adjustment
padj <- mt.rawp2adjp(as.numeric(res.df$pval), proc = "BH")
padj <- padj$adjp[order(padj$index),]
res.df$padj <- padj[,2] # this is the benjamini hochberg adjusted pval
write.csv(res.df, file=here("allele_spec","results", paste0("ase.",vcf_region,".",celltype,".",resample,".output.table")))







