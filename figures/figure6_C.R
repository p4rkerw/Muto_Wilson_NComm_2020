library(Seurat)
library(dplyr)
library(Matrix)
library(biomaRt)
rnaAggr <- readRDS("cellranger_rna_prep/rnaAggr_control.rds")
load("comparison_mouse_datasets/IRIall.RData") #R object from Kirita et.al PNAS_2020
human_gene <- rownames(rnaAggr)
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

converted_genelist = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = human_gene , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
#saveRDS(converted_genelist,"~/Desktop/InterSpecies/converted_genelist.rds")

converted_genelist <- converted_genelist %>% distinct(MGI.symbol,.keep_all=TRUE)

mousePT <- subset(IRI.all,idents = c("PTS1","PTS2","PTS3","NewPT1","NewPT2"))
mousePT_count <- mousePT@assays[["RNA"]]@counts

humanPT <- subset(rnaAggr,idents = c("PT","PT_VCAM1"))
humanPT_count <- as.data.frame(humanPT@assays[["RNA"]]@counts)
humanPT_count$HGNC.symbol <- rownames(humanPT_count)
humanPT_count <-  merge(humanPT_count,converted_genelist,by = "HGNC.symbol")
rownames(humanPT_count) <- humanPT_count$MGI.symbol
humanPT_count <- humanPT_count[,2:5486]
humanPT_count <- as.matrix(humanPT_count)
humanPT_count <- Matrix(humanPT_count, sparse = TRUE)

humanPT <- CreateSeuratObject(counts = humanPT_count)
humanPT <- AddMetaData(humanPT,rnaAggr@meta.data)
humanPT <- NormalizeData(humanPT, assay = 'RNA')
humanPT <- FindVariableFeatures(humanPT, selection.method = "vst", nfeatures = 2000)
humanPT <- ScaleData(humanPT)
humanPT <- RunPCA(humanPT)

mousePT <- CreateSeuratObject(counts = mousePT_count)
mousePT <- AddMetaData(mousePT,IRI.all@meta.data)
mousePT <- NormalizeData(mousePT, assay = 'RNA')
mousePT <- FindVariableFeatures(mousePT, selection.method = "vst", nfeatures = 2000)
mousePT <- ScaleData(mousePT)
mousePT <- RunPCA(mousePT)

PT.anchors <- FindTransferAnchors(
  reference = mousePT, query = humanPT, dims = 1:30)
#saveRDS("Interspecies/PT.anchors.rds")

####################### High resolution ################################

PT <- subset(rnaAggr,idents = c("PT","PT_VCAM1"))
PT <- FindNeighbors(PT, dims = 1:20, verbose = TRUE, reduction = "harmony")
PT <- FindClusters(PT, verbose = TRUE, resolution = 0.1, reduction = "harmony")
PT <- RunUMAP(PT, dims = 1:20, verbose = TRUE, reduction = "harmony")
fig5b_1 <- DimPlot(PT,group.by = "celltype",reduction = "umap",cols = c("#8DCFF7","#F3CB41")) #620x471
predictions <- TransferData(anchorset = PT.anchors, refdata = mousePT$name, dims = 1:30)
PT <- AddMetaData(PT,predictions)
Idents(PT) <- "celltype"
PT1 <- subset(PT,idents = "PT")
PT2 <- subset(PT,idents = "PT_VCAM1")
test1 <- as.data.frame(table(PT1@meta.data[["predicted.id"]]))
colnames(test1) <- c("celltype","PT")
test2 <- as.data.frame(table(PT2@meta.data[["predicted.id"]]))
colnames(test2) <- c("celltype","PT_VCAM1")
test <- merge(test1,test2,by="celltype",all=T)
test[is.na(test)] <- 0
test <-test[c(3,4,5,1,2),]
rownames(test) <- test$celltype
test <- test[,c(2,3)]
fig5b_3 <- pheatmap::pheatmap(test,scale = "column",cluster_rows = F,cluster_cols = F,border_color = "black") #385x628

Idents(PT) <- "predicted.id"
levels(PT) <- c("PTS1","PTS2","PTS3","NewPT2","NewPT1")
new.cluster.ids <- c("PT","PT","PT","FR_PT","PT_Injured")
names(new.cluster.ids) <- levels(PT)
PT <- RenameIdents(PT, new.cluster.ids)

fig5b_2 <- DimPlot(PT,reduction = "umap",cols = c("#8DCFF7","#F3CB41","#F69171")) #620x471