# this script will eliminate doublets from an aggregated snRNA object prior to preprocessing
library(Seurat) # 3.02
library(ggplot2)
library(harmony) # 1.0
library(Rcpp)
library(DoubletFinder)
library(openxlsx)
library(here)
set.seed(1234)

outs <- "cellranger_rna_aggr_control/outs/"

# load aggregated snRNAseq data from cellranger aggregate matrix and create seurat objects
counts <- Read10X(here(outs,"filtered_feature_bc_matrix"))
metadata <- read.csv(here(outs,"aggregation.csv"))

# load aggregated snRNAseq data from cellranger aggregate matrix and create seurat objects
rnaAggr <- CreateSeuratObject(counts = counts, min.cells = 10, min.features = 500, 
                              project = "RNA")

# extract GEM groups from individual barcodes using string split and the suffix integer
# use the GEM groups to assign sample origin (control vs. diabetes) from the aggregation.csv metadata
# files were aggregated control 1, 2, 3, diabetes 1, 2, 3 are correspond to GEM groups 1-6
gemgroup <- sapply(strsplit(rownames(rnaAggr@meta.data), split="-"), "[[", 2) 
current.gemgroups <- seq(1, length(levels(metadata$library_id)))
orig.ident <- levels(metadata$library_id) 
sampleID <- plyr::mapvalues(gemgroup, from = current.gemgroups, to = orig.ident)
rnaAggr <- AddMetaData(object=rnaAggr, metadata=data.frame(orig.ident=sampleID, row.names=rownames(rnaAggr@meta.data)))
rnaAggr <- PercentageFeatureSet(rnaAggr, pattern = "^MT-", col.name = "percent.mt")
rnaAggr <- PercentageFeatureSet(rnaAggr, pattern = "^RPL", col.name = "percent.rpl")
rnaAggr <- PercentageFeatureSet(rnaAggr, pattern = "^RPS", col.name = "percent.rps")

VlnPlot(object = rnaAggr, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "orig.ident",pt.size=0)
VlnPlot(object = rnaAggr, features = c("percent.rps", "percent.rpl"), ncol = 2, group.by = "orig.ident",pt.size=0)
# Removing low quality cells with high mitochondrial RNAs and ribosomal protein RNAs (known to be stable and thus enriched in cells compared to nuclei)

Idents(rnaAggr) <- "orig.ident"

#Dividing by the original ident
Cont1 <- subset(rnaAggr, idents = "Cont1")
Cont2 <- subset(rnaAggr, idents = "Cont2")
Cont3 <- subset(rnaAggr, idents = "Cont3")
Cont4 <- subset(rnaAggr, idents = "Cont4")
Cont5 <- subset(rnaAggr, idents = "Cont5")

# Doublet removal with the assumption that doublets represent 6% of cells.

# Cont1
## Pre-process Seurat object (standard) --------------------------------------------------------------------------------------

VlnPlot(object = Cont1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rps", "percent.rpl"), ncol = 3, group.by = "orig.ident",pt.size=0.1)
# 6464 cells before filtration and 5854 cells after filtration (90.6%)
Cont1 <- subset(Cont1, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & nCount_RNA < 20000 & percent.mt < 2 & percent.rps < 0.8 & percent.rpl < 0.8)
Cont1 <- NormalizeData(Cont1)
Cont1 <- ScaleData(Cont1)
Cont1 <- FindVariableFeatures(Cont1, selection.method = "vst", nfeatures = 2000)
Cont1 <- RunPCA(Cont1)
ElbowPlot(Cont1)
Cont1 <- FindNeighbors(Cont1, dims = 1:20)
Cont1 <- RunUMAP(Cont1, dims = 1:20)
DimPlot(Cont1)

## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list_kidney <- paramSweep_v3(Cont1, PCs = 1:20, sct = F)
sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = FALSE)
bcmvn_kidney <- find.pK(sweep.stats_kidney)

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
Cont1 <- doubletFinder_v3(Cont1, PCs = 1:20, pN = 0.25, pK = 0.02, nExp = round(0.05*length(Cont1@active.ident)), reuse.pANN = FALSE, sct = F)
DimPlot(Cont1,group.by = "DF.classifications_0.25_0.02_293")
FeaturePlot(Cont1,features = "pANN_0.25_0.02_293")

## Doublet Information ----------------------------------------------------------------
doubletdata <- Cont1@meta.data[["DF.classifications_0.25_0.02_293"]]
names(doubletdata) <- rownames(Cont1@meta.data)
doubletdata <- as.data.frame(doubletdata)
Cont1 <- AddMetaData(Cont1,doubletdata)

# Cont2
## Pre-process Seurat object (standard) --------------------------------------------------------------------------------------
VlnPlot(object = Cont2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rps", "percent.rpl"), ncol = 3, group.by = "orig.ident",pt.size=0.1)
# 3600 cells before filtration and 3270 cells after filtration (90.8%)
Cont2 <- subset(Cont2, subset = nFeature_RNA > 500 & nFeature_RNA < 3000 & nCount_RNA < 8000 & percent.mt < 1.5 & percent.rps < 0.8 & percent.rpl < 0.8)
Cont2 <- NormalizeData(Cont2)
Cont2 <- ScaleData(Cont2)
Cont2 <- FindVariableFeatures(Cont2, selection.method = "vst", nfeatures = 2000)
Cont2 <- RunPCA(Cont2)
ElbowPlot(Cont2)
Cont2 <- FindNeighbors(Cont2, dims = 1:17)
Cont2 <- RunUMAP(Cont2, dims = 1:17)
DimPlot(Cont2)

## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list_kidney <- paramSweep_v3(Cont2, PCs = 1:17, sct = F)
sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = FALSE)
bcmvn_kidney <- find.pK(sweep.stats_kidney)
ggplot(bcmvn_kidney, aes(x=pK, y=BCmetric)) + geom_bar(stat = "identity")

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
Cont2 <- doubletFinder_v3(Cont2, PCs = 1:17, pN = 0.25, pK = 0.1, nExp = round(0.05*length(Cont2@active.ident)), reuse.pANN = FALSE, sct = F)
DimPlot(Cont2,group.by = "DF.classifications_0.25_0.1_164")

## Doublet Information ----------------------------------------------------------------
doubletdata <- Cont2@meta.data[["DF.classifications_0.25_0.1_164"]]
names(doubletdata) <- rownames(Cont2@meta.data)
doubletdata <- as.data.frame(doubletdata)
Cont2 <- AddMetaData(Cont2,doubletdata)

# Cont3
## Pre-process Seurat object (standard) --------------------------------------------------------------------------------------

VlnPlot(object = Cont3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rps", "percent.rpl"), ncol = 3, group.by = "orig.ident",pt.size=0.1)
# 6093 cells before filtration and 5764 cells after filtration (94.6%)
Cont3 <- subset(Cont3, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & nCount_RNA < 10000 & percent.mt < 1 & percent.rps < 0.8 & percent.rpl < 0.8)
Cont3 <- NormalizeData(Cont3)
Cont3 <- ScaleData(Cont3)
Cont3 <- FindVariableFeatures(Cont3, selection.method = "vst", nfeatures = 2000)
Cont3 <- RunPCA(Cont3)
ElbowPlot(Cont3)
Cont3 <- FindNeighbors(Cont3, dims = 1:20)
Cont3 <- RunUMAP(Cont3, dims = 1:20)
DimPlot(Cont3)

## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list_kidney <- paramSweep_v3(Cont3, PCs = 1:20, sct = F)
sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = FALSE)
bcmvn_kidney <- find.pK(sweep.stats_kidney)

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
Cont3 <- doubletFinder_v3(Cont3, PCs = 1:20, pN = 0.25, pK = 0.02, nExp = round(0.05*length(Cont3@active.ident)), reuse.pANN = FALSE, sct = F)
DimPlot(Cont3,group.by = "DF.classifications_0.25_0.02_288")
FeaturePlot(Cont3,features = "pANN_0.25_0.02_288")

## Doublet Information ----------------------------------------------------------------
doubletdata <- Cont3@meta.data[["0.25_0.02_288"]]
names(doubletdata) <- rownames(Cont3@meta.data)
doubletdata <- as.data.frame(doubletdata)
Cont3 <- AddMetaData(Cont3,doubletdata)

# Cont4
## Pre-process Seurat object (standard) --------------------------------------------------------------------------------------
VlnPlot(object = Cont4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rps", "percent.rpl"), ncol = 3, group.by = "orig.ident",pt.size=0.1)
# 4048 cells before filtration and 3827 cells after filtration (94.5%)
Cont4 <- subset(Cont4, subset = nFeature_RNA > 500 & nFeature_RNA < 3600 & nCount_RNA < 10000 & percent.mt < 1 & percent.rps < 0.8 & percent.rpl < 0.8)
Cont4 <- NormalizeData(Cont4)
Cont4 <- ScaleData(Cont4)
Cont4 <- FindVariableFeatures(Cont4, selection.method = "vst", nfeatures = 2000)
Cont4 <- RunPCA(Cont4)
ElbowPlot(Cont4)
Cont4 <- FindNeighbors(Cont4, dims = 1:20)
Cont4 <- RunUMAP(Cont4, dims = 1:20)
DimPlot(Cont4)

## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list_kidney <- paramSweep_v3(Cont4, PCs = 1:20, sct = F)
sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = FALSE)
bcmvn_kidney <- find.pK(sweep.stats_kidney)

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
Cont4 <- doubletFinder_v3(Cont4, PCs = 1:20, pN = 0.25, pK = 0.44, nExp = round(0.05*length(Cont4@active.ident)), reuse.pANN = FALSE, sct = F)
DimPlot(Cont4,group.by = "DF.classifications_0.25_0.44_191")

## Doublet Information ----------------------------------------------------------------
doubletdata <- Cont4@meta.data[["DF.classifications_0.25_0.44_191"]]
names(doubletdata) <- rownames(Cont4@meta.data)
doubletdata <- as.data.frame(doubletdata)
Cont4 <- AddMetaData(Cont4,doubletdata)

# Cont5
## Pre-process Seurat object (standard) --------------------------------------------------------------------------------------
VlnPlot(object = Cont5, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rps", "percent.rpl"), ncol = 3, group.by = "orig.ident",pt.size=0.1)
# 4080 cells before filtration and 3933 cells after filtration (95.3%)
Cont5 <- subset(Cont5, subset = nFeature_RNA > 500 & nFeature_RNA < 3000 & nCount_RNA < 7000 & percent.mt < 0.5 & percent.rps < 0.8 & percent.rpl < 0.8)
Cont5 <- NormalizeData(Cont5)
Cont5 <- ScaleData(Cont5)
Cont5 <- FindVariableFeatures(Cont5, selection.method = "vst", nfeatures = 2000)
Cont5 <- RunPCA(Cont5)
ElbowPlot(Cont5)
Cont5 <- FindNeighbors(Cont5, dims = 1:16)
Cont5 <- RunUMAP(Cont5, dims = 1:16)
DimPlot(Cont5)

## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list_kidney <- paramSweep_v3(Cont5, PCs = 1:16, sct = F)
sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = FALSE)
bcmvn_kidney <- find.pK(sweep.stats_kidney)

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
Cont5 <- doubletFinder_v3(Cont5, PCs = 1:16, pN = 0.25, pK = 0.04, nExp = round(0.05*length(Cont5@active.ident)), reuse.pANN = FALSE, sct = F)
DimPlot(Cont5,group.by = "DF.classifications_0.25_0.04_194")

## Doublet Information ----------------------------------------------------------------
doubletdata <- Cont5@meta.data[["DF.classifications_0.25_0.04_194"]]
names(doubletdata) <- rownames(Cont5@meta.data)
doubletdata <- as.data.frame(doubletdata)
Cont5 <- AddMetaData(Cont5,doubletdata)


# Add the doublet data to the original aggregated rna object

doubletdata_1 <- Cont1@meta.data[["doubletdata"]]
doubletdata_2 <- Cont2@meta.data[["doubletdata"]]
doubletdata_3 <- Cont3@meta.data[["doubletdata"]]
doubletdata_4 <- Cont4@meta.data[["doubletdata"]]
doubletdata_5 <- Cont5@meta.data[["doubletdata"]]

names(doubletdata_1) <- rownames(Cont1@meta.data)
names(doubletdata_2) <- rownames(Cont2@meta.data)
names(doubletdata_3) <- rownames(Cont3@meta.data)
names(doubletdata_4) <- rownames(Cont4@meta.data)
names(doubletdata_5) <- rownames(Cont5@meta.data)

doubletdata_1 <- as.data.frame(doubletdata_1)
doubletdata_2 <- as.data.frame(doubletdata_2)
doubletdata_3 <- as.data.frame(doubletdata_3)
doubletdata_4 <- as.data.frame(doubletdata_4)
doubletdata_5 <- as.data.frame(doubletdata_5)

colnames(doubletdata_1) <- "doublet"
colnames(doubletdata_2) <- "doublet"
colnames(doubletdata_3) <- "doublet"
colnames(doubletdata_4) <- "doublet"
colnames(doubletdata_5) <- "doublet"

doubletdata <- rbind(doubletdata_1,doubletdata_2,doubletdata_3,doubletdata_4,doubletdata_5)
dir.create("rna_doublets",showWarnings = FALSE)
saveRDS(doubletdata,here("rna_doublets","rnaAggr_doublets.rds") 

rnaAggr <- AddMetaData(rnaAggr,doubletdata)
Idents(rnaAggr) <- "doublet"
rnaAggr <- subset(rnaAggr,idents = "Singlet")
saveRDS(rnaAggr,here("cellranger_rna_prep","rnaAggr_nodoublets.rds")