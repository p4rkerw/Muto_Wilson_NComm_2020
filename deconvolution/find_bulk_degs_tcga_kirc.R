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



# read in the metadata for the samples and filter for non-tumor KIRC samples
file_meta <- fread("comparison_human_datasets/tcga_kirc/gdc_sample_sheet.2020-05-08.tsv")
file_meta <- dplyr::filter(file_meta, `Sample Type` == "Solid Tissue Normal")
solid_tissue_normal_file <- file_meta$`File ID`
caseID <- file_meta$`Case ID`

# get the demographics for the patients
sample_meta <- fread("comparison_human_datasets/tcga_kirc/sample.tsv")
sample_meta <- sample_meta[sample_meta$case_submitter_id %in% caseID, ]

clinical_meta <- fread("comparison_human_datasets/tcga_kirc/clinical.tsv")
clinical_meta <- clinical_meta[clinical_meta$case_submitter_id %in% caseID, ]
clinical_meta <- dplyr::distinct(clinical_meta, case_id, .keep_all = TRUE)

colnames(file_meta)[1] <- "file_id"
colnames(file_meta)[2] <- "file_name"
colnames(file_meta)[6] <- "case_submitter_id"
metadata <- merge(file_meta, clinical_meta, by="case_submitter_id") %>%
  dplyr::select(file_id, file_name, case_submitter_id, age_at_diagnosis) %>%
  dplyr::mutate(age = as.integer(as.numeric(as.character(age_at_diagnosis))/365.25))

# unzip the files
# lapply(file_paths, function(x) gunzip(paste0("comparison_human_datasets/tcga_kirc_nontumor/",x)))
# grab all kirc file names and paths
files <- list.files("comparison_human_datasets/tcga_kirc", 
                         recursive = TRUE, pattern = "htseq.counts") %>%
  as.data.frame()
colnames(files)[1] <- "file_path"

files$fileID <- str_split(files$file_path, pattern = "/", simplify=TRUE)[,1]
files$file_name <- str_split(files$file_path, pattern = "/", simplify=TRUE)[,2]
files <- files[files$fileID %in% solid_tissue_normal_file, ]
file_paths <- paste0("comparison_human_datasets/tcga_kirc/", files$file_path)

# convert the ensg to hgnc
edb <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
txdb <- getBM(mart = edb, attributes=c("ensembl_gene_id","hgnc_symbol")) # collapse ensembl transcripts to gene id

# read in files
dge.ls <- lapply(file_paths, function(x, TXdb=txdb) {
  df <- fread(x) %>%
    dplyr::rename(gene = V1, counts = V2) %>%
    dplyr::filter(grepl("ENSG", gene)) %>%
    dplyr::arrange(desc(counts, gene)) %>%
    dplyr::mutate(ensg = str_split(gene, pattern = "[.]", simplify = TRUE)[,1]) %>%
    dplyr::select(ensg, counts)
  df <- merge(df, TXdb, by.x="ensg", by.y="ensembl_gene_id") 
  df <- df %>% dplyr::arrange(desc(counts, hgnc_symbol)) %>%
    dplyr::distinct(hgnc_symbol, .keep_all=TRUE) %>%
    dplyr::rename(gene=hgnc_symbol) %>%
    dplyr::select(gene, counts) 
  return(df)
})
names(dge.ls) <- caseID

# generate coldata
coldata <- data.frame(age = metadata$age)
coldata <- dplyr::mutate(coldata, condition = ifelse(age >= 60, "Over_60y","Under_60y")) 
rownames(coldata) <- caseID

# generate a counts matrix
counts <- dge.ls %>% reduce(inner_join, by="gene")
counts <- column_to_rownames(counts, var="gene")
colnames(counts) <- caseID

# import data into deseq2 for rna seq analysis and level the conditions
dds <- DESeqDataSetFromMatrix(counts,
                    colData = coldata,
                    design = ~ condition)

# perform DESeq data analysis and generate vst counts
dds <- DESeq(dds)
vsd <- assay(vst(dds, blind=FALSE))

# convert dds to a dgelist object and then an expressionset to prepare for deconvolution 
# import the vst dds counts for deconvolution
library(DEFormats)
dge <- as.DGEList(dds)
pheno <- as(dge$samples, "AnnotatedDataFrame") # coerces a data.frame to an AnnotatedDataFrame.
eset <- ExpressionSet(assayData = vsd, phenoData = pheno)
saveRDS(eset, file = "comparison_human_datasets/tcga_kirc/eset.tcga.bulkRna.rds")




