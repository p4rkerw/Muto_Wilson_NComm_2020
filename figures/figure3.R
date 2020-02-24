library(Seurat) # 3.0.2
library(Signac) #version 0.1.3
library(BSgenome.Hsapiens.UCSC.hg38)
library(JASPAR2018)
library(TFBSTools)
library(chromVAR)
library(EnsDb.Hsapiens.v86)
library(BuenColors)
set.seed(1234)
library(openxlsx)

#snRNA-seq data
rnaAggr <- readRDS("cellranger_rna_prep/rnaAggr_control.rds")

#Use the object with all version of jasper motifs
sub_atac <- readRDS("cellranger_atac_prep/atacAggr_sub97_control.rds")

marker <- FindAllMarkers(sub_atac, only.pos = T,assay = "chromvar")
marker_set <- tfmarker$gene[!duplicated(tfmarker$gene)]
aver_chromvar <- AverageExpression(sub_atac,assays = "chromvar",features = marker_set)
fig3a <- pheatmap::pheatmap(aver_chromvar[["chromvar"]],scale = "row",
                            cluster_cols=F,cluster_rows = F,
                            color = jdb_palette("brewer_yes"),
                            show_rownames=F)

DefaultAssay(sub_atac) <- "chromvar"
fig3b_1 <- FeaturePlot(sub_atac,features = "HNF4A",cols =jdb_palette("brewer_yes"))
fig3b_2 <- FeaturePlot(sub_atac,features = "TFAP2B",cols =jdb_palette("brewer_yes"))

DefaultAssay(sub_atac) <- "RNA"
fig3b_3 <- FeaturePlot(sub_atac,features = "HNF4A",cols =jdb_palette("Zissou"))
fig3b_4 <- FeaturePlot(sub_atac,features = "TFAP2B",cols =jdb_palette("Zissou"))

#transcriptome data
fig3b_5 <- FeaturePlot(rnaAggr,features = "HNF4A",order=T)
fig3b_6 <- FeaturePlot(rnaAggr,features = "TFAP2B",order=T)



