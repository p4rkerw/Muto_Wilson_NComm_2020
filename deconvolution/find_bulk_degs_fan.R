# load the gene count files from the diabetic nephropathy paper with tximport
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

sample_dir <- list.dirs("comparison_human_datasets/fan_PMID31578193/quants", recursive = FALSE)
files <- file.path(sample_dir, "quant.sf")

# create a metadata df for sample condition
metadata <- str_split(basename(sample_dir), pattern = "_", simplify = TRUE) %>%
  as.data.frame() %>%
  dplyr::select(c(V1,V2,V3)) %>%
  dplyr::mutate(condition = str_sub(V3, end = 1))
colnames(metadata) <- c("Run","GEO_Accession","Sample_ID","condition")

# aggregate transcripts at the gene level using tximport
# import an ensembl transcript database note that some of the transcripts in the reference
# will not be present in the ensembl db because they may represent alternative isoforms
# this has to be adjusted with ignoreTxVersion = TRUE
# files are encoded with "A" for advanced "B" for early DN and "N" for normal
edb <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
txdb <- getBM(mart = edb, attributes=c("ensembl_transcript_id","hgnc_symbol")) # collapse ensembl transcripts to gene id
txi <- tximport(files, type = "salmon", tx2gene = txdb, ignoreTxVersion = TRUE) #ignore version puts all transcripts together

# import data into deseq2 for rna seq analysis and level the conditions
dds <- DESeqDataSetFromTximport(txi,
                                   colData = metadata,
                                   design = ~ condition)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds$condition <- relevel(dds$condition, ref = "N")
levels(dds$condition) <- c("control","advanced","early")

# perform analysis
dds <- DESeq(dds)
res_adv <- results(dds, contrast=c("condition","advanced","control"))
res_early <- results(dds, contrast=c("condition","early","control"))

# perform log-fold shrinkage
resLFC <- lfcShrink(dds, coef="condition_advanced_vs_control", type="apeglm")
resLFC
resOrdered <- res[order(res$pvalue),]
summary(res_early)
sum(res_early$padj < 0.1, na.rm=TRUE)

# visualize results
plotMA(res, ylim=c(-5,5))
plotMA(resLFC, ylim=c(-2,2))

# write out the results
df <- as.data.frame(res_early) %>%
  rownames_to_column(var = "gene") %>%
  dplyr::filter(padj < 0.05)
write.csv(df, file = "bulk_rnaseq/deg_early_bulk.csv", row.names = TRUE)

df <- as.data.frame(res_adv) %>%
  rownames_to_column(var = "gene") %>%
  dplyr::filter(padj < 0.05)
write.csv(df, file = "bulk_rnaseq/deg_adv_bulk.csv", row.names = TRUE)

# plot various genees
plotCounts(dds, gene = "PCK1", intgroup = "condition")

# convert dds to a dgelist object and then an expressionset to prepare for deconvolution 
library(DEFormats)
dge <- as.DGEList(dds)
pheno <- as(dge$samples, "AnnotatedDataFrame") # coerces a data.frame to an AnnotatedDataFrame.
eset <- ExpressionSet(assayData = dge$counts, phenoData = pheno)
saveRDS(eset, file = "bulk_rnaseq/eset.bulkRna.rds")
