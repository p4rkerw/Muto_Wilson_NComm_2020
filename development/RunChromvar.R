# this script will preprocess aggregated snATACseq data from 5 healthy control kidney cortex samples
# counted and aggregated by cellranger-atac v1.2.0 without library normalization

library(Seurat) # 3.0.2
library(Signac) #version 0.2.1
library(BSgenome.Hsapiens.UCSC.hg38)
library(JASPAR2018)
library(TFBSTools)
library(chromVAR)
library(here)
set.seed(1234)
library(openxlsx)

sub_atac <- readRDS("cellranger_atac_prep/atacAggr_sub97_control.rds")

pwm <- getMatrixSet(
  x = JASPAR2018,
  opts = list(species = 9606, all_versions = FALSE)
)

# Scan the DNA sequence of each peak for the presence of each motif
motif.matrix <- CreateMotifMatrix(
  features = StringToGRanges(rownames(sub_atac), sep = c(":", "-")),
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
sub_atac[['peaks']] <- AddMotifObject(
  object = sub_atac[['peaks']],
  motif.object = motif
)

sub_atac <- RegionStats(
  object = sub_atac,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  sep = c(":", "-")
)

# chromVAR can use multiple cores
# library(BiocParallel)
# register(MulticoreParam(8)) # Use 8 cores

# compute motif activities using chromvar
sub_atac <- RunChromVAR(
  object = sub_atac,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  motif.matrix = motif.matrix
)

# Feature name conversion from motifID to TF gene name
x <- NULL                        
for (i in rownames(sub_atac@assays[["chromvar"]]@data)) { 
x <- c(x,pwm@listData[[i]]@name)
}
rownames(sub_atac@assays[["chromvar"]]@data) <- x

saveRDS(sub_atac, file = here("cellranger_atac_prep", "atacAggr_sub97_control.rds"))

