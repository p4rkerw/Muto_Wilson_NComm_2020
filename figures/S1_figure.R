library(Seurat)
library(ggplot2)

rnaAggr <- readRDS("cellranger_rna_prep/rnaAggr_control.rds")
endo <- subset(rnaAggr, idents = c("ENDO"))
endo <- FindNeighbors(endo, dims = 1:20, verbose = TRUE)
endo <- FindClusters(endo, verbose = TRUE, resolution = 0.4)
endo <- RunUMAP(endo, dims = 1:20, verbose = TRUE)

# regroup the idents
endo <- RenameIdents(endo, 
                     '0'='GEC',
                     '1'='AVR',
                     '2'='AVR',
                     '3'='AEA_DVR',
                     '4'='AVR')

p1 <- DimPlot(endo, label=TRUE, repel=TRUE) + 
  NoLegend() +
  ggtitle("A")

# endothelial subcluster markers obtained from Lake et al. PMID: 31249312
markers.to.plot <- c("HECW2","NRG3","FLT1","EMCN", #GEC
                     "DNASE1L3","PECAM1", #AVR
                     "VIM","AQP1" #AEA, DVR
                     )

p2 <- DotPlot(endo, features=rev(markers.to.plot)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab(element_blank()) +
  xlab(element_blank()) +
  ggtitle("B")
CombinePlots(list(p1,p2))
