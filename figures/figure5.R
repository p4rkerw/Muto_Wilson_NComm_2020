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

Fig5a <- DimPlot(tal, reduction = "umap", pt.size = 1) + NoLegend() #610x575

levels(tal) <- rev(new.cluster.ids)
features <- c("SLC12A1","UMOD","CLDN10-AS1","CALCR","CLDN10","CLDN16","CASR","WNK1","JAG1","S100A2","TMPRSS4","S100A6","GADD45B","CRYAB","AKR1B1")
Fig5b <- DotPlot(tal, features = rev(features)) + RotatedAxis()#792x575

Fig5c <- FeaturePlot(tal,features = c("CLDN10","CLDN16","CALCR","CASR"),order=T)  #680x575

#==================ATAC-seq=========================================
library(Signac)
library(GenomeInfoDb)
library(ggplot2)
library(EnsDb.Hsapiens.v86)

sub_atac <- readRDS("cellranger_atac_prep/atacAggr_sub97_control_allT.rds")
tal_atac <- subset(sub_atac, ident = "TAL")
DefaultAssay(tal_atac) <- "peaks"

#tal_atac  <- RunHarmony(tal_atac ,"orig.ident",assay.use = "peaks")
tal_atac   <- RunUMAP(object = tal_atac , reduction = 'harmony', dims = 1:33)
tal_atac   <- FindNeighbors(object = tal_atac , reduction = 'harmony', dims = 1:33)
tal_atac   <- FindClusters(object = tal_atac , verbose = FALSE, resolution = 0.1)

#Annotation of TAL cells
new.cluster.ids <- c("TAL1","TAL2","ATL")
names(new.cluster.ids) <- levels(tal_atac)
tal_atac <- RenameIdents(tal_atac, new.cluster.ids)
tal_atac@meta.data$subtype <- tal_atac@active.ident

fig5e <- DimPlot(tal_atac, pt.size = 1)+NoLegend()  #610x575

DefaultAssay(tal_atac) <- "RNA"
fig5f <- FeaturePlot(tal_atac,features = c("CLDN10","CLDN16","CALCR","CASR"),order=T,cols =jdb_palette("Zissou"))

levels(tal_atac) <- rev(levels(tal_atac))
features <- c("SLC12A1","UMOD","CLDN10-AS1","CALCR","CLDN10","CLDN16","CASR","WNK1","JAG1","S100A2","TMPRSS4","S100A6","GADD45B","CRYAB","AKR1B1")
fig.5g <- DotPlot(tal_atac, features = rev(features)) + RotatedAxis()#792x575

#Motif enrichment analysis on chromvar between TAL1 vs TAL2
DefaultAssay(tal_atac) <- "chromvar"
enriched.motifs_activities <- FindMarkers(tal_atac,ident.1 = "TAL1",ident.2 = "TAL2")

library(ggpubr)
#Visualization
enriched.motifs_activities$motif <- rownames(enriched.motifs_activities)
enriched.motifs_activities$logp <- -log10(enriched.motifs_activities$p_val_adj)
fig.5h <- ggbarplot(head(enriched.motifs_activities), x = "motif", y = "logp",
                    orientation = "horizontal",order=rev(enriched.motifs_activities$motif),
                    fill = "motif", palette = c("#fb6f6f","#fb6f6f","#fb6f6f","#fb6f6f","#5ad40c","#fb6f6f")) #551x575

#Motif enrichment analysis on the dac between TAL1 vs TAL2
DefaultAssay(tal_atac) <- 'peaks'

#Background should be the accessible regions detected at least 2.5% of total TAL population 
background.use <- rownames(tal_atac)[(1-rowSums(tal_atac@assays[["peaks"]]@data==0)/length(colnames(tal_atac)) > 0.025)]
dac <- FindMarkers(tal_atac, ident.1 = "TAL1",ident.2 = "TAL2", test.use = 'LR', latent.vars = "nCount_peaks")

enriched.motifs_TAL1 <- FindMotifs(object = tal_atac, features = rownames(dac[dac$avg_logFC > 0, ]),background = background.use)
enriched.motifs_TAL2 <- FindMotifs(object = tal_atac, features = rownames(dac[dac$avg_logFC < 0, ]),background = background.use)

#Visualization
enriched.motifs_TAL1$logp <- -log10(enriched.motifs_TAL1$pvalue)
enriched.motifs_TAL2$logp <- -log10(enriched.motifs_TAL2$pvalue)

fig.5i <- ggbarplot(head(enriched.motifs_TAL1), x = "motif.name", y = "logp",
                    orientation = "horizontal",order=rev(head(enriched.motifs_TAL1)$motif.name),
                    fill = "motif.name", palette = c("#5ad40c","#5ad40c","#5ad40c","#5ad40c","#5ad40c","#5ad40c")) #551x575

fig.5j <- ggbarplot(head(enriched.motifs_TAL2), x = "motif.name", y = "logp",
                    orientation = "horizontal",order=rev(head(enriched.motifs_TAL2)$motif.name),
                    fill = "motif.name", palette = c("#fb6f6f","#fb6f6f","#fb6f6f","#fb6f6f","#fb6f6f","#fb6f6f")) #551x575


