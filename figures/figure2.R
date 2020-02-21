library(Signac) #version 0.1.5
library(Seurat) #version 3.0.2
library(GenomeInfoDb)
library(harmony) #version 1.0
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(here)
library(tibble)
library(dplyr)
library(openxlsx)
library(readxl)
set.seed(1234)

sub_atac <- readRDS("cellranger_atac_prep/atacAggr_sub97_control.rds")
#dar.celltype.control.xlsx is put in the cellranger_atac_prep folder
darfile <- "cellranger_atac_prep/dar.celltype.control.xlsx"
x <- NULL                        
for (i in 1:14) {x <-  rbind(x, read.xlsx(darfile,sheet = i))
}                             
dar_total <- x
dar_pos <- subset(dar_total,dar_total$avg_logFC > 0)
peaklist <- dar_pos[,1]
peaklist <- peaklist[!duplicated(peaklist)]
DefaultAssay(sub_atac) <- "peaks"
Idents(sub_atac) <- "celltype"
dac_aver <- AverageExpression(sub_atac,features = peaklist,assays = "peaks")
fig2a <- pheatmap::pheatmap(dac_aver[["peaks"]],scale='row',cluster_rows=F,cluster_cols=F)

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

x <- NULL                     
for (i in 1:14) {x <-  c(x, list(read.xlsx(file,sheet = i)))
}                   
darlist <- x
darlist <- c(darlist,list(dar_total))

for (i in 1:15) {darlist[[i]] <-  StringToGRanges(regions=darlist[[i]][,1], sep = c(":", "-"))
}                   

names(darlist) <- celltype <- c("PCT","PST","PT_KIM1","PEC","TAL",
                                "DCT","CNT","PC","ICA","ICB",
                                "PODO","ENDO","MES_FIB","LEUK","TOTAL")

peakAnnoList <- lapply(darlist, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE)

fig2c <- plotAnnoPie(peakAnnoList[["TOTAL"]]) #total DAR in the dataset
fig2d <- plotAnnoBar(peakAnnoList[1:14]) #celltype-specific analysis
fig2e_1 <- plotDistToTSS(peakAnnoList[1:14])

promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
CELLtag <- getTagMatrix(darlist[["TOTAL"]], windows=promoter)
fig2e_2 <- tagHeatmap(CELLtag, xlim=c(-3000, 3000), color="#3587B6")
fig2e_3 <- plotAvgProf(CELLtag, xlim=c(-3000, 3000), conf = 0.95, resample = 1000)




