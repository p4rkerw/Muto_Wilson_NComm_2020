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
library(matrixStats)

#snRNA-seq data
rnaAggr <- readRDS("cellranger_rna_prep/rnaAggr_control.rds")

#Use the object with all version of jasper motifs
sub_atac <- readRDS("cellranger_atac_prep/chromVar.atacAggr_sub97_control.rds")

# find cell type specific TF and make a unique list
tf_celltype.df <- FindAllMarkers(sub_atac, only.pos = T, assay = "chromvar")
unique_tf.list <- unique(tf_celltype.df$gene)
aver_chromvar <- AverageExpression(sub_atac, assays = "chromvar", features = unique_tf.list) %>%
  as.data.frame()

# sort the heatmap prior to plotting find column index for maxium value for each TF activity
aver_chromvar <- aver_chromvar[do.call(order, c(aver_chromvar, list(decreasing=TRUE))),]
aver_chromvar$max <- max.col(aver_chromvar)
aver_chromvar <- aver_chromvar[order(aver_chromvar$max), ]
aver_chromvar <- dplyr::select(aver_chromvar, -max)
colnames(aver_chromvar) <- idents

# visualize results
fig3a <- pheatmap::pheatmap(aver_chromvar,scale = "row",
                            cluster_cols=F,cluster_rows = F,
                            color = jdb_palette("brewer_yes"),
                            show_rownames=F)

# HNF4A does not appear to be in the pwm for JASPAR2018
DefaultAssay(sub_atac) <- "chromvar"
fig3b_1 <- FeaturePlot(sub_atac,features = "HNF4A",cols =jdb_palette("brewer_yes"))
fig3b_2 <- FeaturePlot(sub_atac,features = "TFAP2B",cols =jdb_palette("brewer_yes"))

# HNF1B and TFAP2B
fig3b_1 <- FeaturePlot(sub_atac,features = "MA0153.2",cols =jdb_palette("brewer_yes")) + ggtitle("HNF1B")
fig3b_2 <- FeaturePlot(sub_atac,features = "MA0811.1",cols =jdb_palette("brewer_yes")) + ggtitle("TFAP2B")

DefaultAssay(sub_atac) <- "RNA"
fig3b_3 <- FeaturePlot(sub_atac,features = "HNF4A",cols =jdb_palette("Zissou"))
fig3b_4 <- FeaturePlot(sub_atac,features = "TFAP2B",cols =jdb_palette("Zissou"))

#transcriptome data
fig3b_5 <- FeaturePlot(rnaAggr,features = "HNF4A",order=T)
fig3b_6 <- FeaturePlot(rnaAggr,features = "TFAP2B",order=T)



