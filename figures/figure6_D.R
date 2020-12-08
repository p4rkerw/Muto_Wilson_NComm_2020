# devtools::install_github('shenorrlab/bseqsc')
library(bseqsc)
library(here)
library(BisqueRNA) # for Seurat to ExpressionSet
library(Seurat)
library(dplyr)
library(ggplot2)
library(reshape)
library(ggpubr)
library(BuenColors)

rnaAggr <- readRDS("cellranger_rna_prep/rnaAggr_control.rds")

Idents(rnaAggr) <- "celltype"
new.cluster.ids <- c("PT","PT_VCAM1","PEC","TAL","DCT_CNT",
                     "DCT_CNT","DCT_CNT","PC","ICA","ICB",
                     "PODO","ENDO","MES_FIB","MES_FIB","LEUK")
names(new.cluster.ids) <- levels(rnaAggr)
rnaAggr <- RenameIdents(rnaAggr, new.cluster.ids)
rnaAggr@meta.data$celltype <- rnaAggr@active.ident

scrna_eset <- SeuratToExpressionSet(rnaAggr, delimiter="-", position=2, version="v3")
scrna_eset@phenoData$celltype <- rnaAggr$celltype

# read in bulk eset
bulk_eset <- readRDS("comparison_human_datasets/tcga_kirc/eset.tcga.bulkRna.rds")

# add estimated celltype proportions as phenoData
res <- BisqueRNA::ReferenceBasedDecomposition(bulk_eset, scrna_eset, markers=NULL, use.overlap=F)
pData(bulk_eset) <- cbind(pData(bulk_eset), t(res[["bulk.props"]]))


# metadata located here
toplot <- bulk_eset@phenoData@data %>%
  dplyr::select(group, PT_VCAM1)

toplot$TCGA <- "TCGA"

p2 <- ggplot(toplot, aes(y=PT_VCAM1, fill=TCGA)) + 
  geom_boxplot() + 
  xlab("Patient Age Range") +
  ylab("Proportion of PT_VCAM1 Cells") +
  ggtitle("Proportion of PT_VCAM1 Cells by Patient Age Range TCGA-KIRC")


fig5d <- ggboxplot(toplot, x = "TCGA", y = "PT_VCAM1",
                   add = "jitter",ylim = c(0, 0.2),
fill = "TCGA",palette = "#A3CDE5")+NoLegend() 


