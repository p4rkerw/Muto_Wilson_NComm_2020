# load the gene count files from TCGA with tximport
library(readr)
library(dplyr)
library(stringr)
library(data.table)
library(DESeq2)
library(biomaRt)
library(tibble)
library(GEOquery)
library(purrr)
# library(EnsDb.Hsapiens.v86)

# read in the metadata for the samples and filter for non-tumor samples
file_meta <- fread("comparison_human_datasets/tcga_kidney_nontumor/gdc_sample_sheet.2020-05-12.tsv")
file_meta <- dplyr::filter(file_meta, `Sample Type` == "Solid Tissue Normal")
solid_tissue_normal_file <- file_meta$`File ID`
caseID <- file_meta$`Case ID`

# get the demographics for the patients
sample_meta <- fread("comparison_human_datasets/tcga_kidney_nontumor/sample.tsv")
sample_meta <- sample_meta[sample_meta$case_submitter_id %in% caseID, ]

clinical_meta <- fread("comparison_human_datasets/tcga_kidney_nontumor/clinical.tsv")
clinical_meta <- clinical_meta[clinical_meta$case_submitter_id %in% caseID, ]
clinical_meta <- dplyr::distinct(clinical_meta, case_id, .keep_all = TRUE)

colnames(file_meta)[1] <- "file_id"
colnames(file_meta)[2] <- "file_name"
colnames(file_meta)[6] <- "case_submitter_id"
metadata <- merge(file_meta, clinical_meta, by="case_submitter_id") %>%
  dplyr::select(file_id, file_name, case_submitter_id, age_at_diagnosis) %>%
  dplyr::mutate(age = as.integer(as.numeric(as.character(age_at_diagnosis))/365.25))

# unzip the files
# lapply(file_paths, function(x) gunzip(paste0("comparison_human_datasets/tcga_kidney_nontumor/",x)))

# grab all kirc file names and paths
files <- list.files("comparison_human_datasets/tcga_kidney_nontumor/", 
                         recursive = TRUE, pattern = "htseq.counts") %>%
  as.data.frame()
colnames(files)[1] <- "file_path"

files$fileID <- str_split(files$file_path, pattern = "/", simplify=TRUE)[,1]
files$file_name <- str_split(files$file_path, pattern = "/", simplify=TRUE)[,2]
files <- files[files$fileID %in% solid_tissue_normal_file, ]
file_paths <- paste0("comparison_human_datasets/tcga_kidney_nontumor/", files$file_path)

# read in files
dge.ls <- lapply(file_paths, function(x) {
  df <- fread(x) %>%
    dplyr::rename(gene = V1, counts = V2) %>%
    dplyr::filter(grepl("ENSG", gene)) %>%
    dplyr::arrange(desc(counts, gene)) %>%
    dplyr::mutate(ensg = str_split(gene, pattern = "[.]", simplify = TRUE)[,1]) %>%
    dplyr::select(ensg, counts)
  return(df)
})
names(dge.ls) <- caseID

# generate coldata
coldata <- data.frame(age = metadata$age)
coldata <- dplyr::mutate(coldata, condition = ifelse(age >= 60, "Over_60y","Under_60y")) 
rownames(coldata) <- caseID


# generate a counts matrix
counts <- dge.ls %>% reduce(full_join, by="ensg") %>%
  dplyr::arrange(ensg)
counts <- column_to_rownames(counts, var="ensg")
colnames(counts) <- caseID

# import data into deseq2 for rna seq analysis and level the conditions
dds <- DESeqDataSetFromMatrix(counts,
                    colData = coldata,
                    design = ~ condition)

# perform DESeq data analysis and extract normalized counts
dds <- DESeq(dds)
counts <- data.frame(counts(dds, normalized=TRUE)[-1,]) %>% # remove the normalization factor in the first row
  rownames_to_column(var="ensg")

# convert the ensg to hgnc
edb <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
txdb <- getBM(mart = edb, attributes=c("ensembl_gene_id","hgnc_symbol")) # collapse ensembl transcripts to gene id
counts <- merge(counts, txdb, by.x="ensg", by.y="ensembl_gene_id")
counts$mean <- rowMeans(as.matrix(counts[2:129]), na.rm = TRUE) # avg count per gene
counts <- dplyr::arrange(counts, desc(hgnc_symbol, mean)) %>%
  dplyr::distinct(hgnc_symbol, .keep_all=TRUE) %>%
  dplyr::select(-ensg, -mean) %>%
  column_to_rownames(var="hgnc_symbol")

# convert dds to a dgelist object and then an expressionset to prepare for deconvolution 
library(DEFormats)
dge <- as.DGEList(dds)
colnames(counts) <- rownames(dge$samples)
pheno <- as(dge$samples, "AnnotatedDataFrame") # coerces a data.frame to an AnnotatedDataFrame.
eset <- ExpressionSet(assayData = as.matrix(counts), phenoData = pheno)
saveRDS(eset, file = "comparison_human_datasets/tcga_kidney_nontumor/eset.tcga.nontumor.bulkRna.rds")




