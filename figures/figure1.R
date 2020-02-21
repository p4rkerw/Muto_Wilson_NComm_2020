library(Signac) #version 0.1.5
library(Seurat) #version 3.0.2
library(GenomeInfoDb)
library(harmony) #version 1.0
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(here)
library(tibble)
library(dplyr)
set.seed(1234)

#snRNA-seq

rnaAggr <- readRDS("cellranger_rna_prep/rnaAggr_control.rds")
Idents(rnaAggr) <- "celltype"

fig1a <- DimPlot(rnaAggr, reduction = "umap",label = T, repel = T) + NoLegend()

# reorder the idents and save celltype annotations in celltype slot metadata
levels(rnaAggr) <- rev(levels(rnaAggr))
features <- c("SLC34A1","HAVCR1","CFH","SLC12A1","SLC12A3",
              "SLC8A1","AQP2","SLC26A7","SLC26A4","NPHS2",
              "EMCN","PIEZO2","COL1A2","PTPRC")
fig1b <- DotPlot(rnaAggr, features = rev(features), cols = c("lightyellow","royalblue")) +
  RotatedAxis() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


#snATAC-seq
sub_atac <- readRDS("cellranger_atac_prep/atacAggr_sub97_control.rds")
atacAggr <- readRDS("cellranger_atac_prep/atacAggr_control.rds")

Idents(sub_atac) <- "highres.predicted.id"
fig1d_1 <- DimPlot(sub_atac,label = T,repel = T) + NoLegend()
fig1c_5 <- DimPlot(sub_atac, label = T, repel = T,group.by = "highres.predicted.id") + NoLegend()


Idents(sub_atac) <- "celltype"
fig1d_2 <- DimPlot(sub_atac, label = T, repel = T) + NoLegend()

# reorder the idents and save celltype annotations in celltype slot metadata
levels(sub_atac) <- rev(c("PCT","PST","PT_KIM1","PEC","TAL",
                          "DCT","CNT","PC","ICA","ICB",
                          "PODO","ENDO","MES_FIB","LEUK"))
features <- c("SLC34A1","SLC5A2","SLC5A1","HAVCR1","CFH",
              "SLC12A1","SLC12A3","SLC8A1","AQP2","SLC26A7",
              "SLC26A4","NPHS2","EMCN","ACTA2","PTPRC")
DefaultAssay(sub_atac) <- "RNA"
fig1e <- DotPlot(sub_atac, features = rev(features),cols = c("lightyellow","royalblue")) +
  RotatedAxis()


rnaAggr@meta.data[["tech"]] <- "RNA"
atacAggr@meta.data[["tech"]] <- "ATAC"
#hue_pal()(2) 
#"#F8766D" "#00BFC4"
Idents(rnaAggr) <- "tech"
fig1c_1 <- DimPlot(rnaAggr,cols = "#00BFC4",reduction = "umap") + NoLegend()
fig1c_2 <- DimPlot(atacAggr,cols = "#F8766D",reduction = "umap") + NoLegend()

#Use the transfer.anchor for label-transfer in signacAtacAggrPreprocess
rnaAggr <- NormalizeData(rnaAggr,assay = "RNA")
DefaultAssay(rnaAggr) <- "RNA"
transfer.anchors <- readRDS(here("cellranger_atac_prep","transfer_anchors_control.rds"))
genes.use <- transfer.anchors@anchor.features
refdata <- GetAssayData(rnaAggr, assay = "RNA", slot = "data")[genes.use, ]
imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = atacAggr[["pca"]])
atacAggr[["RNA"]] <- imputation
DefaultAssay(atacAggr) <- "RNA"
coembed <- merge(x = rnaAggr, y = atacAggr)
# Finally, we run PCA and UMAP on this combined object, to visualize the co-embedding of both
DefaultAssay(coembed) <- "RNA"
coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE,assay = "RNA")
coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE,assay = "RNA")
coembed <- RunHarmony(coembed, "orig.ident", plot_convergence = TRUE, assay.use = "RNA")
coembed <- RunUMAP(coembed, dims = 1:20, reduction = 'harmony')
coembed$celltype1 <- ifelse(!is.na(coembed$highres.predicted.id), coembed$highres.predicted.id, coembed$celltype)
Idents(coembed) <- "celltype1" # these are the snRNA celltype annotations
levels(coembed) <- c("PT","PT_KIM1","PEC","TAL","DCT1",
                     "DCT2","CNT","PC","ICA","ICB",
                     "PODO","ENDO","MES","FIB","LEUK")

fig1c_3 <- DimPlot(coembed,group.by = "tech")
fig1c_4 <- DimPlot(coembed,label = T,repel = T) + NoLegend()













