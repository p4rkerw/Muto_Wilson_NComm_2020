# load the gene count files from Shavlakadze et al with tximport
library(tximport)
library(readr)
library(dplyr)
library(stringr)
library(data.table)
library(DESeq2)
library(biomaRt)
library(tibble)
# library(EnsDb.Hsapiens.v86)

sample_dir <- list.dirs("comparison_mouse_datasets/shavlakadze_PMID31533046/salmon_rat", recursive = FALSE)
files <- file.path(sample_dir, "quant.sf")

# create a metadata df for sample condition
metadata <- data.frame(basename(sample_dir)) %>%
  dplyr::mutate(condition = str_split(basename.sample_dir., pattern = "_", simplify=TRUE)[,1]) %>%
  dplyr::mutate(replicate = str_split(basename.sample_dir., pattern = "_", simplify=TRUE)[,2]) 

# aggregate transcripts at the gene level using tximport
rat <- useMart("ensembl", dataset = "rnorvegicus_gene_ensembl")
txdb <- getBM(mart = rat, attributes=c("ensembl_transcript_id","rgd_symbol")) # collapse ensembl transcripts to gene id
txi <- tximport(files, type = "salmon", tx2gene = txdb, ignoreTxVersion = TRUE) #ignore version puts all transcripts together

# import data into deseq2 for rna seq analysis and level the conditions
dds <- DESeqDataSetFromTximport(txi,
                                   colData = metadata,
                                   design = ~ condition)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# perform DESeq data analysis and extract normalized counts
dds <- DESeq(dds)
counts <- counts(dds, normalized=TRUE)[-1,] # remove the normalization factor in the first row

# convert rat genes to human hgnc
ratGenes = rownames(counts)
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
lookup_table = getLDS(attributes = c("rgd_symbol"), filters = "rgd_symbol", values = ratGenes , mart = rat,
                      attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)

# extract the counts from genes that are shared between species
colnames(lookup_table)[1] <- "symbol"
counts <- data.frame(assay(dds))[-1,]
counts <- rownames_to_column(counts, var = "symbol")
counts <- merge(x=lookup_table, y=counts, by="symbol") 
counts$rowMean <- rowMeans(as.matrix(counts[3:ncol(counts)]), na.rm = TRUE)
counts <- dplyr::arrange(counts, desc(HGNC.symbol, rowMean)) %>%
  dplyr::distinct(HGNC.symbol, .keep_all=TRUE) %>%
  dplyr::select(-symbol, -rowMean) %>%
  column_to_rownames(var="HGNC.symbol")

# convert dds to a dgelist object and then an expressionset to prepare for deconvolution 
library(DEFormats)
dge <- as.DGEList(dds)
colnames(counts) <- rownames(dge$samples)
pheno <- as(dge$samples, "AnnotatedDataFrame") # coerces a data.frame to an AnnotatedDataFrame.
eset <- ExpressionSet(assayData = as.matrix(counts), phenoData = pheno)
saveRDS(eset, file = "comparison_mouse_datasets/shavlakadze_PMID31533046/eset.shavlakadze.bulkRna.rds")
