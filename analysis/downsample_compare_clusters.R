# this script will downsample the snRNA and snATAC datasets to contain the same number of cells from each donor
# and compare the ability of snRNA and snATAC to resolve different cell types

library(Seurat)
library(here)
library(dplyr)
library(tibble)
library(harmony)
library(Signac)
library(matrixStats)
library(reticulate)
use_python(python = "~/Anaconda3/bin/python", required = TRUE)

rnaAggr <- readRDS(here("cellranger_rna_prep","rnaAggr_control.rds"))
num_rna <- table(rnaAggr@meta.data$orig.ident) %>% as.data.frame()

atacAggr <- readRDS(here("cellranger_atac_prep","atacAggr_sub97_control.rds"))
num_atac <- table(atacAggr@meta.data$orig.ident) %>% as.data.frame()

# the number of atac cells is greater in each donor 
merge <- data.frame(rna=num_rna$Freq, atac=num_atac$Freq, donor = num_atac$Var1)

# downsample the aggregated snATAC object using the corresponding number of snRNA cells per donor
donors <- merge$donor
Idents(atacAggr) <- "orig.ident"
set.seed(1234)
bc_keep.ls <- lapply(donors, function(donor) {
  meta <- rownames_to_column(atacAggr@meta.data, var="barcode") %>%
    dplyr::filter(orig.ident == donor) 
  num_keep <- merge[merge$donor ==  donor,]$rna
  bc_keep <- sample(meta$barcode, size=num_keep, replace=F)
}) 
bc_keep <- unlist(bc_keep.ls)

# subset the aggregated snATAC object
atacSub <- subset(atacAggr, cells = bc_keep)

# process the subset
DefaultAssay(atacSub) <- "peaks"
atacSub <- RunTFIDF(atacSub) #TF-IDF normalization again in the subset data.
atacSub <- FindTopFeatures(atacSub, min.cutoff = 'q1')
atacSub <- RunSVD(
  object = atacSub,
  assay = 'peaks',
  reduction.key = 'pca_', # this is actually an LSI reduction called "pca"
  reduction.name = 'pca')
atacSub <- RunHarmony(atacSub, "orig.ident", plot_convergence = TRUE, assay.use = "peaks")
atacSub <- RunUMAP(object = atacSub, reduction = 'harmony', dims = 1:29, assay.use = "peaks")
atacSub <- FindNeighbors(object = atacSub, reduction = 'harmony', dims = 1:29, assay.use = "peaks")
atacSub <- FindClusters(object = atacSub, verbose = FALSE, reduction = 'harmony', assay.use = "peaks")
DimPlot(atacSub)

num_atac_sub <- table(atacSub@meta.data$orig.ident) %>% as.data.frame()

# the number of atac cells is greater in each donor 
merge <- data.frame(rna=num_rna$Freq, atac=num_atac_sub$Freq, donor = num_atac_sub$Var1)

# compare distribution of cells
atac_dist <- table(atacAggr@meta.data$orig.ident, atacAggr@meta.data$celltype)
total_col = rowSums(atac_dist)
pcts <- round(atac_dist / rowSums(atac_dist), digits=2)

rna_dist <- table(rnaAggr@meta.data$orig.ident, rnaAggr@meta.data$celltype)
total_col = rowSums(rna_dist)
pcts <- round(rna_dist / rowSums(rna_dist), digits=2)



