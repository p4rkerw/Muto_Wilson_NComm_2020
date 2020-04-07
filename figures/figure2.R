library(Signac) #version 0.2.1
library(Seurat) #version 3.0.2
library(GenomeInfoDb)
library(harmony) #version 1.0
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(here)
library(tibble)
library(dplyr)
library(openxlsx)
set.seed(1234)

sub_atac <- readRDS("cellranger_atac_prep/atacAggr_sub97_control.rds")
darfile <- "analysis_control/dar.celltype.control.xlsx"
idents <- getSheetNames(darfile)
list.dar <- lapply(idents, function(x) {
  df <- read.xlsx(darfile, sheet = x, rowNames = T) %>%
    rownames_to_column(var = "coord") %>%
    dplyr::mutate(celltype = x) # annotate each region with its corresponding celltype
  })

# identify all unique cell-type-specific peaks and filter for logfc > 0
all_dar <- bind_rows(list.dar) %>%
  dplyr::filter(avg_logFC > 0) %>%
  dplyr::select("coord") %>%
  dplyr::distinct()

dar_aver <- AverageExpression(sub_atac, features = all_dar$coord, assays = "peaks")
fig2a <- pheatmap::pheatmap(dar_aver[["peaks"]], scale = 'row', cluster_rows = F, cluster_cols = F,
                            show_rownames = FALSE)

sub_atac <- SetFragments(object = sub_atac, file = here(outs_atac, "fragments.tsv.gz"))

fig2b <- CoveragePlot(
  object = sub_atac,
  region = "chr2:169359547-169365451",
  sep = c(":", "-"),
  extend.upstream = 5000,
  extend.downstream = 5000,
  annotation = EnsDb.Hsapiens.v86,
)


library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
library(clusterProfiler)
library(org.Hs.eg.db)

# convert the DAR to GRanges objects to annotate
all_dar.gr <- StringToGRanges(all_dar$coord, sep = c(":","-"))
list.dar.gr <- lapply(seq(list.dar), function(x) {
  df <- list.dar[[x]]
  gr <- StringToGRanges(df$coord, sep = c(":","-"))
  return(gr)
})
names(list.dar.gr) <- idents

# annotate the list of GRanges DAR for each cell type
list.peakAnno <- lapply(list.dar.gr, annotatePeak, TxDb = txdb,
                       tssRegion = c(-3000, 3000), verbose = FALSE)
all.peakAnno <- annotatePeak(all_dar.gr, TxDb = txdb,
                       tssRegion = c(-3000, 3000), verbose = FALSE)

fig2c <- plotAnnoPie(all.peakAnno) #total DAR in the dataset
fig2d <- plotAnnoBar(list.peakAnno) #celltype-specific analysis
fig2e_1 <- plotDistToTSS(list.peakAnno)

promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
CELLtag <- getTagMatrix(all.peakAnno, windows=promoter)
fig2e_2 <- tagHeatmap(CELLtag, xlim=c(-3000, 3000), color="#3587B6")
fig2e_3 <- plotAvgProf(CELLtag, xlim=c(-3000, 3000), conf = 0.95, resample = 1000)




