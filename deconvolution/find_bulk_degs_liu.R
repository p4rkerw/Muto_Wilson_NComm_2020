# load the gene count files from Liu et al PMID28931758 
library(readr)
library(dplyr)
library(stringr)
library(data.table)
library(DESeq2)
library(biomaRt)
library(tibble)
library(GEOquery)
library(purrr)
library(openxlsx)
# library(EnsDb.Hsapiens.v86)

# read in the FPKM counts 
mouseCounts <- read.xlsx("comparison_mouse_datasets/liu_PMID28931758/GSE98622_mouse-iri-master.xlsx") %>%
  dplyr::filter(grepl("ENSM", X1)) %>%
  dplyr::distinct(symbol, .keep_all = TRUE)

# convert mouse genes to human genes
mouseGenes <- mouseCounts$symbol
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
lookup_table = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = mouseGenes , mart = mouse,
                   attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)

# extract the counts from genes that are shared between species
colnames(lookup_table)[1] <- "symbol"
counts <- merge(x=lookup_table, y=mouseCounts, by="symbol") 
counts$mean <- rowMeans(as.matrix(counts[5:ncol(counts)]), na.rm = TRUE)
counts <- dplyr::group_by(counts, HGNC.symbol) %>%
  dplyr::arrange(desc(mean), .by_group = TRUE) %>%
  dplyr::distinct(HGNC.symbol, .keep_all=TRUE) %>%
  dplyr::select(-X1, -type, -symbol, -mean) %>%
  column_to_rownames(var="HGNC.symbol")

# convert to integer counts
i <- seq(colnames(counts))
counts[ , i] <- apply(counts[ , i], 2,            
                    function(x) as.integer(x))

# generate coldata
coldata <- data.frame(sample = colnames(counts))
coldata$condition <- str_split(coldata$sample, pattern = "-", simplify = TRUE)[,1]

# import data into deseq2 for rna seq analysis and level the conditions
dds <- DESeqDataSetFromMatrix(counts,
                    colData = coldata,
                    design = ~ condition)

# do not put through DESeq2 because these are FPKM as opposed to read counts or transcript abundance
# dds <- DESeq(dds)
# counts <- counts(dds, normalized=TRUE)[-1]
# do not normalize the FPKM values

# convert dds to a dgelist object and then an expressionset to prepare for deconvolution 
library(DEFormats)
dge <- as.DGEList(dds)
pheno <- as(dge$samples, "AnnotatedDataFrame") # coerces a data.frame to an AnnotatedDataFrame.
eset <- ExpressionSet(assayData = dge$counts, phenoData = pheno)
saveRDS(eset, file = "comparison_mouse_datasets/liu_PMID28931758/eset.liu.bulkRna.rds")




