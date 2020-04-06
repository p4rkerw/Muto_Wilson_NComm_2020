library(Seurat)

atacAggr <- readRDS("cellranger_atac_prep/atacAggr_control.rds")
hist(atacAggr@meta.data$prediction.score.max,
     main="",
     xlab="Prediction Score",
     col="red",
     xlim=c(0,1))
