# load the gene count files from Cippa et al with tximport
library(tximport)
library(readr)
library(dplyr)
library(stringr)
library(data.table)
library(DESeq2)
library(biomaRt)
library(tibble)
# library(EnsDb.Hsapiens.v86)

# import the quant.sf files output by salmon 
# dir <- system.file("extdata", package="tximportData")
# samples <- read.table(file.path(dir,"samples.txt"), header=TRUE)
# samples$condition <- factor(rep(c("A","B"),each=3))
# rownames(samples) <- samples$run
# samples[,c("pop","center","run","condition")]

sample_dir <- list.dirs("comparison_human_datasets/cippa_PMID30429361/quants", recursive = FALSE)
files <- file.path(sample_dir, "quant.sf")

# create a metadata df for sample condition
metadata <- str_split(basename(sample_dir), pattern = "[.]", simplify = TRUE) %>%
  as.data.frame() %>%
  dplyr::mutate(V2 = gsub(pattern="_baseline", replacement = "", V2)) %>%
  dplyr::mutate(Sample_ID = str_split(V2, pattern="_", simplify=TRUE)[,2]) %>%
  dplyr::mutate(condition = str_split(V2, pattern="_", simplify=TRUE)[,3]) %>%
  dplyr::mutate(GEO_Accession = str_split(V1, pattern="_", simplify=TRUE)[,1]) %>%
  dplyr::select(Sample_ID, condition, GEO_Accession)

# aggregate transcripts at the gene level using tximport
# import an ensembl transcript database note that some of the transcripts in the reference
edb <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
txdb <- getBM(mart = edb, attributes=c("ensembl_transcript_id","hgnc_symbol")) # collapse ensembl transcripts to gene id
txi <- tximport(files, type = "salmon", tx2gene = txdb, ignoreTxVersion = TRUE) #ignore version puts all transcripts together

# import data into deseq2 for rna seq analysis and level the conditions
dds <- DESeqDataSetFromTximport(txi,
                                   colData = metadata,
                                   design = ~ condition)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds$condition <- relevel(dds$condition, ref = "pre")
levels(dds$condition) <- c("pre","post","3months","1year")

# perform DESeq data analysis and extract normalized counts 
dds <- DESeq(dds)
counts <- counts(dds, normalized=TRUE)[-1,] # remove the normalization factor in the first row

# convert dds to a dgelist object and then an expressionset to prepare for deconvolution 
library(DEFormats)
dge <- as.DGEList(dds)
colnames(counts) <- rownames(dge$samples)
pheno <- as(dge$samples, "AnnotatedDataFrame") # coerces a data.frame to an AnnotatedDataFrame.
eset <- ExpressionSet(assayData = as.matrix(counts), phenoData = pheno)
saveRDS(eset, file = "comparison_human_datasets/cippa_PMID30429361/eset.cippa.bulkRna.rds")
