library(Seurat)
library(Signac)
library(ggplot2)

rnaAggr <- readRDS("cellranger_rna_prep/rnaAggr_control.rds")
pt_rna <- subset(rnaAggr, idents = "PT")
p1 <- FeaturePlot(pt_rna, features=c("SLC5A1")) +
  ggtitle("A) SLC5A1 Expression")
p2 <- FeaturePlot(pt_rna, features=c("SLC5A2")) +
  ggtitle("B) SLC5A2 Expression")


atacAggr <- readRDS("cellranger_atac_prep/atacAggr_sub97_control.rds")
pt_atac <- subset(atacAggr, idents = c("PCT","PST"))
p3 <- FeaturePlot(pt_atac, features=c("SLC5A1")) +
  ggtitle("C) SLC5A1 Gene Activity")
p4 <- FeaturePlot(pt_atac, features=c("SLC5A2")) +
  ggtitle("D) SLC5A2 Gene Activity")


CombinePlots(list(p1,p2,p3,p4), ncol=2)
