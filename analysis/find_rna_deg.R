# this script will take an aggregated seurat object and identify differentially expressed genes
#  between celltypes

library(Seurat) # 3.0.2
library(openxlsx)
library(here)

rnaAggr <- readRDS(here("cellranger_rna_prep","rnaAggr_control.rds"))
Idents(rnaAggr) <- "celltype"
DefaultAssay(rnaAggr) <- 'SCT'

# wrapper function for FindMarkers 
GetMarkers <- function(cluster, seurat_aggregate) {
  print(paste0("Finding DEG for: ",cluster))
  rnaAggr <- seurat_aggregate
  deg <- FindMarkers(rnaAggr, 
                     ident.1 = cluster,    
                     min.pct = 0.2) # find all cluster-specific degs
  return(deg)
}

# FindMarkers and write to an xlsx file with default parameters
idents <- levels(rnaAggr@meta.data$celltype)
list.cluster.deg <- lapply(idents, function(x) {GetMarkers(x, seurat_aggregate = rnaAggr)})

dir.create(here("analysis_control"), showWarnings = FALSE)
write.xlsx(list.cluster.deg, file = here("analysis_control","deg.celltype.control.xlsx"), sheetName = idents, rowNames = T)

