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

Idents(rnaAggr) <- "seurat_clusters"
new.cluster.ids <- c("PT","DCT1","TAL","CNT","PT",
                     "TAL","ICA","TAL","PC","ENDO",
                     "PEC","DCT2","PODO","DR-PT","ICB",
                     "ENDO","MES","FIB","ENDO","LEUK")
names(new.cluster.ids) <- levels(rnaAggr)
rnaAggr <- RenameIdents(rnaAggr, new.cluster.ids)
levels(rnaAggr) <- c("PT","DR-PT","PEC","TAL","DCT1",
                     "DCT2","CNT","PC","ICA","ICB",
                     "PODO","ENDO","MES","FIB","LEUK")
rnaAggr@meta.data[["celltype"]] <- rnaAggr@active.assay

fig1a <- DimPlot(rnaAggr,reduction = "umap",label = T,repel = T)+NoLegend()

# reorder the idents and save celltype annotations in celltype slot metadata
levels(rnaAggr) <- rev(c("PT","DR-PT","PEC","TAL","DCT1",
                         "DCT2","CNT","PC","ICA","ICB",
                         "PODO","ENDO","MES","FIB","LEUK"))

features <- c("SLC34A1","HAVCR1","CFH","SLC12A1","SLC12A3","SLC8A1","AQP2","SLC26A7","SLC26A4","NPHS2","EMCN","PIEZO2","COL1A2","PTPRC")
fig1b <- DotPlot(rnaAggr, features = rev(features),cols = c("lightyellow","royalblue")) + RotatedAxis()


#snATAC-seq
######For now this code uses dev ver atacAggr_sub95_control(signacAggrPreprocess_yoshi.R)#######
sub_atac <- readRDS("cellranger_atac_prep/atacAggr_sub95_control_dev.rds")
atacAggr <- readRDS("cellranger_atac_prep/atacAggr_control.rds")

Idents(sub_atac) <- "highres.predicted.id"

current.name<-c("PCT","PT_KIM1","PEC","TAL","DCT1","DCT2","CNT","PC","ICA","ICB","PODO","ENDO","MES","FIB","LEUK")
new.name<-c("PT","DR-PT","PEC","TAL","DCT1","DCT2","CNT","PC","ICA","ICB","PODO","ENDO","MES","FIB","LEUK")
sub_atac@active.ident <- plyr::mapvalues(sub_atac@active.ident,from = current.name,to=new.name)


levels(sub_atac) <- c("PT","DR-PT","PEC","TAL","DCT1",
                     "DCT2","CNT","PC","ICA","ICB",
                     "PODO","ENDO","MES","FIB","LEUK")
fig1d_1 <- DimPlot(sub_atac,label = T,repel = T)+NoLegend()


Idents(sub_atac) <- "seurat_clusters"
new.cluster.ids <- c("PCT","TAL","PST","PCT","DCT",
                     "TAL","TAL","PCT","PC","CNT",
                     "ENDO","ICB","ICA","DR-PT","PEC",
                     "DCT","MES-FIB","ENDO","PODO","LEUK")
names(new.cluster.ids) <- levels(sub_atac)
sub_atac <- RenameIdents(sub_atac, new.cluster.ids)

levels(sub_atac) <- c("PCT","DR-PT","PEC","TAL","DCT","PST",
                          "CNT","PC","ICA","ICB","PODO"
                          ,"ENDO","MES-FIB","LEUK")

sub_atac@meta.data$celltype <- sub_atac@active.ident

fig1d_2 <- DimPlot(sub_atac,label = T,repel = T)+NoLegend()

# reorder the idents and save celltype annotations in celltype slot metadata
levels(sub_atac) <- rev(c("PCT","PST","DR-PT","PEC","TAL","DCT",
                     "CNT","PC","ICA","ICB","PODO"
                     ,"ENDO","MES-FIB","LEUK"))
features <- c("SLC34A1","SLC5A2","SLC5A1","HAVCR1","CFH","SLC12A1","SLC12A3","SLC8A1","AQP2","SLC26A7","SLC26A4","NPHS2","EMCN","ACTA2","PTPRC")
fig1e <- DotPlot(sub_atac, features = rev(features),cols = c("lightyellow","royalblue")) + RotatedAxis() 


rnaAggr@meta.data[["tech"]] <- "RNA"
atacAggr@meta.data[["tech"]] <- "ATAC"
#hue_pal()(2) 
#"#F8766D" "#00BFC4"
Idents(rnaAggr) <- "tech"
Idents(atacAggr) <- "tech"
fig1c_1 <- DimPlot(rnaAggr,cols = "#00BFC4",reduction = "umap")+NoLegend()
fig1c_2 <- DimPlot(atacAggr,cols = "#F8766D",reduction = "umap")+NoLegend()

#Use the transfer.anchor for label-transfer in signacAtacAggrPreprocess
genes.use <- transfer.anchors@anchor.features
refdata <- GetAssayData(rnaAggr, assay = "RNA", slot = "data")[genes.use, ]
imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = atacAggr[["lsi"]])
atacAggr[["SCT"]] <- imputation
coembed <- merge(x = rnaAggr, y = atacAggr)

# Finally, we run PCA and UMAP on this combined object, to visualize the co-embedding of both

coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE,assay = "SCT")
coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE,assay = "SCT")
coembed <- RunHarmony(coembed, "orig.ident", plot_convergence = TRUE, assay.use = "SCT")
coembed <- RunUMAP(coembed, dims = 1:24, reduction = 'harmony')
coembed$celltype1 <- ifelse(!is.na(coembed$celltype), coembed$celltype, coembed$highres.predicted.id)
Idents(coembed) <- "celltype1"
levels(coembed) <- c("PCT","PT_KIM1","PEC","TAL","DCT1",
                     "DCT2","CNT","PC","ICA","ICB",
                     "PODO","ENDO","MES","FIB","LEUK")

fig1c_3 <- p5_merge1b <- DimPlot(coembed, group.by = "tech")
fig1c_4 <- DimPlot(coembed,label = T,repel = T)+NoLegend()

fig1c_5 <- hist(atacAggr@meta.data$prediction.score.max)











