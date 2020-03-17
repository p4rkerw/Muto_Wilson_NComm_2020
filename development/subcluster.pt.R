library(Seurat)
library(Signac)

atacAggr <- readRDS("cellranger_atac_prep/atacAggr_sub97_control.rds")
pt <- subset(atacAggr, idents = c("PCT","PST"))
DefaultAssay(pt) <- "peaks"
pt <- RunTFIDF(pt)
pt <- FindTopFeatures(pt, min.cutoff = 'q0')
pt <- RunSVD(pt, assay = 'peaks', reduction.key = 'LSI_', reduction.name = 'lsi')
pt <- RunUMAP(pt, dims = 1:20, verbose = TRUE, reduction = "lsi")
pt <- FindNeighbors(pt, dims = 1:20, verbose = TRUE)
pt <- FindClusters(pt, verbose = TRUE, resolution = 0.7)
DimPlot(pt, label=TRUE)