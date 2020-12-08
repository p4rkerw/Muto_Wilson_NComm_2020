# this script will take an aggregated seurat object and identify differentially accessible chromatin
# regions between celltypes

library(Seurat) # 3.0.2
library(Signac) #version 0.2.1
library(EnsDb.Hsapiens.v86)
library(openxlsx)
library(here)

atacAggr <- readRDS(here("cellranger_atac_prep","atacAggr_sub97_control.rds"))
Idents(atacAggr) <- "celltype"
DefaultAssay(atacAggr) <- 'peaks'

# wrapper function for FindMarkers 
GetMarkers <- function(cluster, seurat_aggregate) {
  print(paste0("Finding DAR for: ",cluster))
  atacAggr <- seurat_aggregate
  dar <- FindMarkers(atacAggr, 
                     ident.1 = cluster,    
                     test.use = 'LR', 
                     latent.vars = "peak_region_fragments",
                     min.pct = 0.2) # find all cluster-specific dars
  cf <- ClosestFeature(rownames(dar), annotation=EnsDb.Hsapiens.v86, sep=c(':','-'))
  return(cbind(dar, gene=cf$gene_name, distance=cf$distance))
}

# FindMarkers and write to an xlsx file with default parameters
idents <- levels(atacAggr@meta.data$celltype)
list.cluster.dar <- lapply(idents, function(x) {GetMarkers(x, seurat_aggregate = atacAggr)})

dir.create(here("analysis_control"), showWarnings = FALSE)
write.xlsx(list.cluster.dar, file = here("analysis_control","dar.celltype.control.xlsx"), sheetName = idents, rowNames = T)

