library(DoubletFinder)
library(Seurat)
library(here)
library(ggplot2)

# retrieve library ids to process individually and identify doublets
# to exclude from aggregated object
rna_counts <- "cellranger_rna_counts/version_3.1.0"
libs <- dir(here(rna_counts))

kidney.data <- here(rna_counts,libs[1],"outs","filtered_feature_bc_matrix")
counts <- Read10X(kidney.data)

## Pre-process Seurat object (sctransform) -----------------------------------------------------------------------------------
seu_kidney <- CreateSeuratObject(counts, min.cells = 10, min.features = 500, project = "RNA")
seu_kidney <- PercentageFeatureSet(seu_kidney, pattern = "^MT-", col.name = "percent.mt")
seu_kidney <- PercentageFeatureSet(seu_kidney, pattern = "^RPL", col.name = "percent.rpl")
seu_kidney <- PercentageFeatureSet(seu_kidney, pattern = "^RPS", col.name = "percent.rps")
seu_kidney <- subset(seu_kidney, 
                     subset = nFeature_RNA > 500 
                      & nFeature_RNA < 6000 
                      & nCount_RNA < 16000 
                      & percent.mt < 0.8 
                      & percent.rps < 0.4 
                      & percent.rpl < 0.4)
seu_kidney <- SCTransform(seu_kidney)
seu_kidney <- RunPCA(seu_kidney)
seu_kidney <- RunUMAP(seu_kidney, dims = 1:12)
# ElbowPlot(seu_kidney, ndims = 40) # to determin number of dimensions for clustering
seu_kidney <- FindNeighbors(seu_kidney, dims = 1:20, verbose = TRUE, reduction = "pca")
seu_kidney <- FindClusters(seu_kidney, verbose = TRUE, resolution = 0.6, reduction = "pca")

## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list_kidney <- paramSweep_v3(seu_kidney, PCs = 1:12, sct = TRUE)
sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = FALSE)
bcmvn_kidney <- find.pK(sweep.stats_kidney)

## pK Identification (ground-truth) ------------------------------------------------------------------------------------------
# sweep.res.list_kidney <- paramSweep_v3(seu_kidney, PCs = 1:12, sct = TRUE)
# gt.calls <- seu_kidney@meta.data[rownames(sweep.res.list_kidney[[1]]), "GT"]
# sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = TRUE, GT.calls = gt.calls)
# bcmvn_kidney <- find.pK(sweep.stats_kidney)

# ? Add in label transfer for automatic annotation of celltypes

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- seu_kidney@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.075*length(seu_kidney@active.ident))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
seu_kidney <- doubletFinder_v3(seu_kidney, PCs = 1:12, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)

Idents(seu_kidney) <- "DF.classifications_0.25_0.09_346"
DimPlot(seu_kidney, reduction="umap", pt.size=0.5)

# identify what celltype the doublets are
p1 <- DimPlot(seu_kidney, reduction = "umap", assay = "SCT", label = TRUE) + ggtitle("Original Clustering")

celltype.markers <- c("CUBN","LRP2","HAVCR1","SLC5A1","SLC5A2", # PT and PT-KIM1+ markers
                      "CFH", # PEC
                      "SLC12A1", # TAL NKCC2
                      "SLC12A3","TRPM6", # DCT1 and DCT2 NCC
                      "SCNN1G","TRPV5", # DCT2/CNT ENaC
                      "CALB1", # CNT
                      "AQP2", # PC
                      "ATP6V0D2", # ICA and ICB
                      "SLC4A1", # ICA
                      "SLC26A4", # ICB
                      "NPHS1","NPHS2", # PODO
                      "PECAM1","FLT1", # ENDO
                      "ITGA8","PDGFRB", # MES
                      "PTPRC") # WBC
p2 <- DotPlot(seu_kidney, features = celltype.markers) + theme(axis.text.x = element_text(angle = 90, hjust = 1))











