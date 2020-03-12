library(Seurat)
library(harmony)
library(ggplot2)
library(dplyr)
library(Matrix)
library(BuenColors)
library(here)
set.seed(1234)

#snRNA-seq data
rnaAggr <- readRDS("cellranger_rna_prep/rnaAggr_control.rds")

#Sub-clustering
tal <- subset(rnaAggr, ident = "TAL")
tal<- FindNeighbors(tal, dims = 1:25, verbose = TRUE, reduction = "harmony")
tal <- FindClusters(tal, verbose = TRUE, resolution = 0.2, reduction = "harmony")
tal <- RunUMAP(tal, dims = 1:25, verbose = TRUE, reduction = "harmony")
DimPlot(tal, reduction = "umap", label = TRUE, pt.size = 1) + NoLegend()

#VlnPlot(tal,features = c("nCount_RNA","nFeature_RNA","percent.mt"),pt.size = 0)
#cluster1 is low-quality cells (low "nCount_RNA","nFeature_RNA" and high "percent.mt") and removed for the sub-clustering downstream analysis of TAL
tal <- subset(tal, idents = 1, invert = T)

#Re-clustering the TAL without low-quality cells
tal<- FindNeighbors(tal, dims = 1:20, verbose = TRUE, reduction = "harmony")
tal <- FindClusters(tal, verbose = TRUE, resolution = 0.2, reduction = "harmony")
tal <- RunUMAP(tal, dims = 1:20, verbose = TRUE, reduction = "harmony")

#Annotation of TAL cells
new.cluster.ids <- c("TAL1","TAL2","ATL")
names(new.cluster.ids) <- levels(tal)
tal <- RenameIdents(tal, new.cluster.ids)
tal@meta.data$subtype <- tal@active.ident

Fig5a <- DimPlot(tal, reduction = "umap", label = TRUE, pt.size = 1) + NoLegend()

levels(tal) <- rev(new.cluster.ids)
features <- c("SLC12A1","UMOD","PHACTR1","CLDN10-AS1","CLDN10","CLDN16","JAG1","WNK1","CASR","S100A2","TMPRSS4","S100A6","GADD45B","CRYAB","AKR1B1")
Fig5b <- DotPlot(tal, features = rev(features)) + RotatedAxis()

Fig5c <- FeaturePlot(tal,features = c("CLDN10","JAG1","CLDN16","NOTCH2"),order=T)

#==================ATAC-seq=========================================
library(Signac)
library(GenomeInfoDb)
library(ggplot2)
library(EnsDb.Hsapiens.v86)

sub_atac <- readRDS("cellranger_atac_prep/sub_atac_sub97_control_allT.rds")
tal_atac <- subset(sub_atac, ident = "TAL")
DefaultAssay(tal_atac) <- "peaks"

#tal_atac <- RunTFIDF(tal_atac)
#tal_atac <- FindTopFeatures(tal_atac, min.cutoff = 'q1')
#tal_atac <- RunSVD(
  #object = tal_atac,
  #assay = 'peaks',
  #reduction.key = 'pca_',
  #reduction.name = 'pca')

#tal_atac  <- RunHarmony(tal_atac ,"orig.ident",assay.use = "peaks")
tal_atac   <- RunUMAP(object = tal_atac , reduction = 'harmony', dims = 1:33)
tal_atac   <- FindNeighbors(object = tal_atac , reduction = 'harmony', dims = 1:33)
tal_atac   <- FindClusters(object = tal_atac , verbose = FALSE, resolution = 0.1)

#Annotation of TAL cells
new.cluster.ids <- c("TAL1","TAL2","ATL")
names(new.cluster.ids) <- levels(tal_atac)
tal_atac <- RenameIdents(tal_atac, new.cluster.ids)
tal_atac@meta.data$subtype <- tal_atac@active.ident

fig5e <- DimPlot(tal_atac,label = T)+NoLegend()

fig5f <- FeaturePlot(tal_atac,features = c("CLDN10","JAG1","CLDN16","NOTCH2"),order=T,cols =jdb_palette("Zissou"))

DimPlot(tal_atac,label=T)
DefaultAssay(tal_atac) <- "RNA"

library(BSgenome.Hsapiens.UCSC.hg38)
library(JASPAR2018)
library(TFBSTools)
library(chromVAR)
set.seed(1234)

DefaultAssay(tal_atac) <- "peaks"

pwm <- getMatrixSet(
  x = JASPAR2018,
  opts = list(species = 9606, all_versions = T)
)

# Scan the DNA sequence of each peak for the presence of each motif
motif.matrix <- CreateMotifMatrix(
  features = StringToGRanges(rownames(tal_atac), sep = c(":", "-")),
  pwm = pwm,
  genome = 'BSgenome.Hsapiens.UCSC.hg38',
  sep = c(":", "-")
)

# Create a new Mofif object to store the results
motif <- CreateMotifObject(
  data = motif.matrix,
  pwm = pwm
)

# Add the Motif object to the assay
tal_atac[['peaks']] <- AddMotifObject(
  object = tal_atac[['peaks']],
  motif.object = motif
)

tal_atac <- RegionStats(
  object = tal_atac,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  sep = c(":", "-")
)

# compute motif activities using chromvar
tal_atac <- RunChromVAR(
  object = tal_atac,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  motif.matrix = motif.matrix
)

DefaultAssay(tal_atac) <- 'chromvar'

# Feature name conversion from motifID to TF gene name
x <- NULL                        
for (i in rownames(tal_atac@assays[["chromvar"]]@data)) { 
  x <- c(x,pwm@listData[[i]]@name)
}
rownames(tal_atac@assays[["chromvar"]]@data) <- x


motifidlist <- rownames(tal_atac@assays[["chromvar"]]@data)

fig6h <- VlnPlot(sub_atac,features = c("RBPJ"),pt.size = 0,ncol=1)+NoLegend()

