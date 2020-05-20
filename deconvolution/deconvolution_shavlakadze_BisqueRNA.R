library(here)
library(BisqueRNA)
library(Seurat)
library(dplyr)
library(ggplot2)
library(reshape)

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
bulk_eset <- readRDS("comparison_mouse_datasets/shavlakadze_PMID31533046/eset.shavlakadze.bulkRna.rds")

res <- BisqueRNA::ReferenceBasedDecomposition(bulk_eset, scrna_eset, markers=NULL, use.overlap=F)
pData(bulk_eset) <- cbind(pData(bulk_eset), t(res[["bulk.props"]]))

# metadata located here
bulk_eset@phenoData@data$group <- factor(bulk_eset@phenoData@data$group, levels=c("6m","9m","12m","18m","21m","24m","27m"))

toplot <- bulk_eset@phenoData@data %>%
  dplyr::select(group, PT_VCAM1)
p1 <- ggplot(toplot, aes(x=group, y=PT_VCAM1, fill=group)) + 
  geom_boxplot() + 
  geom_point(position=position_dodge(width=0.75),aes(group=group)) +
  xlab("Rat Age") +
  ylab("Proportion of PT_VCAM1 Cells") +
  ggtitle("Proportion of PT_VCAM1 Cells by Rat Age")

toplot <- bulk_eset@phenoData@data %>%
  dplyr::select(group, PT)
p2 <- ggplot(toplot, aes(x=group, y=PT, fill=group)) + 
  geom_boxplot() + 
  geom_point(position=position_dodge(width=0.75),aes(group=group)) +
  xlab("Rat Age") +
  ylab("Proportion of PT Cells") +
  ggtitle("Proportion of PT Cells by Rat Age")

CombinePlots(list(p1,p2))

toplot <- bulk_eset@phenoData@data
toplot <- toplot[7:ncol(toplot)]
toplot <- melt(toplot)
p3 <- ggplot(toplot, aes(x=variable, y=value, fill=variable)) + 
  geom_boxplot() + 
  ylab("Proportion of Cells") +
  ggtitle("Proportion of Rat Cells in all Groups") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


toplot <- bulk_eset@phenoData@data %>%
  dplyr::select(group, c(PT,PT_VCAM1))
toplot$rate <- toplot$PT_VCAM1/toplot$PT
p4 <- ggplot(toplot, aes(x=group, y=rate, fill=group)) + 
  geom_boxplot() +
  ggtitle("PT_VCAM1/Normal_PT ratio (Rat data)")

