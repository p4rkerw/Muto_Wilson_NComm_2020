
library(Seurat) # 3.0.2
library(Signac) #version 0.2.1
library(BSgenome.Hsapiens.UCSC.hg38)
library(JASPAR2018)
library(TFBSTools)
library(chromVAR)
set.seed(1234)
library(openxlsx)
library(here)

atacAggr <- readRDS(here("cellranger_atac_prep","atacAggr_sub97_control.rds"))
Idents(atacAggr) <- "celltype"
DefaultAssay(atacAggr) <- 'peaks'

# Get a list of motif position weight matrices from the JASPAR database
pwm <- getMatrixSet(
  x = JASPAR2018,
  opts = list(species = 9606, all_versions = FALSE)
)

# Scan the DNA sequence of each peak for the presence of each motif
motif.matrix <- CreateMotifMatrix(
  features = StringToGRanges(rownames(atacAggr), sep = c(":", "-")),
  pwm = pwm,
  genome = 'BSgenome.Hsapiens.UCSC.hg38',
  sep = c(":", "-")
)

# Create a new Mofif object to store the results
motif <- CreateMotifObject(
  data = motif.matrix,
  pwm = pwm
)

# Add the Motif object to the assay
atacAggr[['peaks']] <- AddMotifObject(
  object = atacAggr[['peaks']],
  motif.object = motif
)

atacAggr <- RegionStats(
  object = atacAggr,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  sep = c(":", "-")
)

# chromVAR can use multiple cores
# library(BiocParallel)
# register(MulticoreParam(8)) # Use 8 cores

# compute motif activities using chromvar
atacAggr <- RunChromVAR(
  object = atacAggr,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  motif.matrix = motif.matrix
)

saveRDS(atacAggr, file = "cellranger_atac_prep/chromVar.atacAggr_sub97_control.rds")

GetMotifs <- function(cluster, seurat_aggregate) {
  print(paste0("Finding motifs for: ",cluster))
  atacAggr <- seurat_aggregate
  DefaultAssay(atacAggr) <- 'peaks'
  dac <- FindMarkers(atacAggr,
                     ident.1 = cluster,
                     test.use = 'LR',
                     latent.vars = "nCount_peaks") # find all cluster-specific degs
  enriched.motifs <- FindMotifs(object = atacAggr, features = rownames(dac[dac$p_val < 0.05, ]))
  return(enriched.motifs)
}

GetChromvarActivities <- function(cluster, seurat_aggregate, motif) {
  print(paste0("Finding chromVAR activities for: ",cluster))
  atacAggr <- seurat_aggregate
  DefaultAssay(atacAggr) <- 'chromvar'
  dam <- FindMarkers(atacAggr, 
                     ident.1 = cluster,
                     test.use = 'LR',
                     latent.vars = "nCount_peaks",
                     logfc.threshold = 0) # find all cluster-specific degs
  motifLookup <- rownames(dam)
  motifNames <- sapply(motifLookup, function(x) motif@motif.names[[x]])
  return(cbind(dam, gene = motifNames))
}

# FindMarkers and write to an xlsx file with default parameters
Idents(atacAggr) <- "celltype"
idents <- levels(atacAggr)
list.cluster.dac <- lapply(idents, function(x) GetMotifs(x, seurat_aggregate = atacAggr))
write.xlsx(list.cluster.dac, file = "analysis_control/motifs.celltype.control.xlsx", sheetName = idents, rowNames = T)

list.cluster.dam <- lapply(idents, function(x) GetChromvarActivities(x, seurat_aggregate = atacAggr, motif = motif))
write.xlsx(list.cluster.dam, file = "analysis_control/chromVAR.celltype.control.xlsx", sheetName = idents, rowNames = T)
