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

atacAggr <- readRDS("cellranger_atac_prep/atacAggr_control.rds")
rnaAggr <- readRDS("cellranger_rna_prep/rnaAggr_control.rds")


# identify anchors to transfer cell labels from snRNAseq to snATACseq "SCT" gene activity scores
transfer.anchors2 <- FindTransferAnchors(
  reference = rnaAggr,
  query = atacAggr,
  reference.assay = 'SCT',
  query.assay = 'RNA',
  reduction = 'cca',
)

lowres.predicted.labels_SCT <- TransferData(
  anchorset = transfer.anchors2,
  refdata = rnaAggr$lowres.celltype,
  weight.reduction = atacAggr[["lsi"]]
)

colnames(lowres.predicted.labels_SCT) <- paste("lowres",colnames(lowres.predicted.labels_SCT), sep = ".")
colnames(lowres.predicted.labels_SCT) <- paste("ver2",colnames(lowres.predicted.labels_SCT), sep = "_")
atacAggr <- AddMetaData(atacAggr, metadata = lowres.predicted.labels_SCT)

table(atacAggr$prediction.score.max > 0.95)
#FALSE  TRUE 
#5579 26533

table(atacAggr$ver2_lowres.prediction.score.max > 0.95)
#FALSE  TRUE 
#5395 26717 


