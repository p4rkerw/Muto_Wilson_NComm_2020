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

fig1a <- DimPlot(rnaAggr, reduction = "umap") + NoLegend() #720x580
figs3a <- DimPlot(rnaAggr,reduction = "umap",split.by = "orig.ident") #500x1500

# reorder the idents and save celltype annotations in celltype slot metadata
levels(rnaAggr) <- rev(levels(rnaAggr))
features <- c("SLC34A1","LRP2","HAVCR1","CFH","SLC12A1","SLC12A3",
              "SLC8A1","AQP2","SLC26A7","SLC26A4","NPHS2",
              "EMCN","PIEZO2","COL1A2","PTPRC")
fig1b <- DotPlot(rnaAggr, features = rev(features), cols = c("lightyellow","royalblue")) +
  RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) #820x580

Idents(rnaAggr) <- "seurat_clusters"
figs1a <- DimPlot(rnaAggr, reduction = "umap",label = T) + NoLegend() #720x580
levels(rnaAggr) <- rev(c(0,4,13,10,2,5,7,1,11,3,8,6,14,12,9,15,18,16,17,19))
features <- c("SLC34A1","LRP2","SLC5A1","SLC5A2","HAVCR1","CFH","SLC12A1","SLC12A3",
              "SLC8A1","AQP2","SLC26A7","SLC26A4","NPHS2",
              "FLT1","EMCN","PECAM1","PIEZO2","COL1A2","PTPRC")

figs1b <- DotPlot(rnaAggr, features = rev(features), cols = c("lightyellow","royalblue")) +
  RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) #820x600


#snATAC-seq
sub_atac <- readRDS("cellranger_atac_prep/atacAggr_sub97_control.rds")
atacAggr <- readRDS("cellranger_atac_prep/atacAggr_control.rds")

Idents(sub_atac) <- "highres.predicted.id"
fig1d_1 <- DimPlot(sub_atac,label = T,repel = T) + NoLegend()
fig1c_5 <- DimPlot(sub_atac, label = T, repel = T,group.by = "highres.predicted.id") + NoLegend()
figs3b <- DimPlot(sub_atac,reduction = "umap",split.by = "orig.ident") #500x1500


Idents(sub_atac) <- "celltype"
fig1d_2 <- DimPlot(sub_atac) + NoLegend() #720x580

# reorder the idents and save celltype annotations in celltype slot metadata
levels(sub_atac) <- rev(levels(sub_atac))
features <- c("SLC34A1","LRP2","SLC5A2","SLC5A1","HAVCR1","CFH",
              "SLC12A1","SLC12A3","SLC8A1","AQP2","SLC26A7",
              "SLC26A4","NPHS2","EMCN","ACTA2","PTPRC")
DefaultAssay(sub_atac) <- "RNA"
fig1e <- DotPlot(sub_atac, features = rev(features),cols = c("lightyellow","royalblue")) +
  RotatedAxis() #820x580


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
levels(coembed) <- c("PT","PT_VCAM1","PEC","TAL","DCT1",
                     "DCT2","CNT","PC","ICA","ICB",
                     "PODO","ENDO","MES","FIB","LEUK")

fig1c_3 <- DimPlot(coembed,group.by = "tech")
fig1c_4 <- DimPlot(coembed,label = T,repel = T) + NoLegend()

Idents(sub_atac) <- "seurat_clusters"
figs2b <- DimPlot(sub_atac, reduction = "umap",label = T) + NoLegend() #720x580
levels(sub_atac) <- rev(c(3,0,2,7,11,15,4,1,6,14,16,5,9,8,13,12,19,10,17,18))
features <- c("SLC34A1","LRP2","SLC5A1","SLC5A2","HAVCR1","CFH","SLC12A1","SLC12A3",
              "SLC8A1","AQP2","SLC26A7","SLC26A4","NPHS2",
              "FLT1","EMCN","PECAM1","PIEZO2","COL1A2","PTPRC")

figs2c <- DotPlot(sub_atac, features = rev(features), cols = c("lightyellow","royalblue")) +
  RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) #820x600












