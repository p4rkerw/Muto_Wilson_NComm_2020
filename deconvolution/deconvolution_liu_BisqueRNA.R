# devtools::install_github('shenorrlab/bseqsc')
library(bseqsc)
library(here)
library(BisqueRNA) # for Seurat to ExpressionSet
library(Seurat)
library(dplyr)
library(ggplot2)
library(reshape)


# load the bulk eset
bulk_eset <- readRDS("comparison_mouse_datasets/liu_PMID28931758/eset.liu.bulkRna.rds")
pData(bulk_eset)$group <- gsub(x=pData(bulk_eset)$group, pattern="SHAM", replacement='SHAM_')
pData(bulk_eset)$group <- gsub(x=pData(bulk_eset)$group, pattern="NORM", replacement='NORM_')
pData(bulk_eset)$group <- gsub(x=pData(bulk_eset)$group, pattern="IRI", replacement='IRI_')
pData(bulk_eset)$group <- gsub(x=pData(bulk_eset)$group, pattern="6mN", replacement='6m') # mislabeled


# read in aggregated snRNA dataset
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

# build reference basis matrix of expression profiles for individual celltypes
# and estimate cell proportions
res <- BisqueRNA::ReferenceBasedDecomposition(bulk_eset, scrna_eset, markers=NULL, use.overlap=F)
pData(bulk_eset) <- cbind(pData(bulk_eset), t(res[["bulk.props"]]))

# metadata located here
bulk_eset@phenoData@data$group <- factor(bulk_eset@phenoData@data$group, 
                                            levels=c("NORM_3m","NORM_9m","NORM_15m","SHAM_4h","SHAM_24h","SHAM_12m",
                                                     "IRI_2h","IRI_4h","IRI_24h","IRI_48h","IRI_72h",
                                                     "IRI_7d","IRI_14d","IRI_28d","IRI_6m","IRI_12m"))


# determine if there is a statistically significant difference between mean PT_VCAM1
# proportion across time points using one-way ANOVA
toplot <- bulk_eset@phenoData@data
res.aov <- aov(PT_VCAM1 ~ group, data=toplot)
summary(res.aov) # p=0.0051

toplot <- bulk_eset@phenoData@data
res.aov <- aov(PT ~ group, data=toplot)
summary(res.aov) # p=0.0051

sham <- pData(bulk_eset)[grepl(x=pData(bulk_eset)$group, pattern="SHAM"),]
iri <- pData(bulk_eset)[grepl(x=pData(bulk_eset)$group, pattern="IRI"),]


toplot <- bulk_eset@phenoData@data %>%
  dplyr::select(group, PT_VCAM1) 
toplot <- toplot[grepl("SHAM|IRI", toplot$group),]
p2 <- ggplot(toplot, aes(x=group, y=PT_VCAM1, fill=group)) + 
  geom_boxplot() + 
  xlab("Time Post-Perfusion") +
  ylab("PT_VCAM1") +
  ggtitle("Estimated Proportion of PT_VCAM1 \nMurine IRI (Liu 2017) (n=49)") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  NoLegend()


toplot <- bulk_eset@phenoData@data %>%
  dplyr::select(group, PT)
toplot <- toplot[grepl("SHAM|IRI", toplot$group),]
p1 <- ggplot(toplot, aes(x=group, y=PT, fill=group)) + 
  geom_boxplot() + 
  xlab("Time Post-Perfusion") +
  ylab("PT") +
  ggtitle("Estimated Proportion of PT \nMurine IRI (Liu 2017) (n=49)") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  NoLegend()

CombinePlots(list(p1,p2))

toplot <- bulk_eset@phenoData@data
toplot <- toplot[5:ncol(toplot)]
toplot <- melt(toplot)
p3 <- ggplot(toplot, aes(x=variable, y=value, fill=variable)) + 
  geom_boxplot() + 
  ylab("Proportion of Cells") +
  ggtitle("Proportion of Cells in Liu et al.")

toplot <- bulk_eset@phenoData@data %>%
  dplyr::filter(group == "IRI_6m")
toplot <- toplot[5:ncol(toplot)]
toplot <- melt(toplot)
p3 <- ggplot(toplot, aes(x=variable, y=value, fill=variable)) + 
  geom_boxplot() + 
  ylab("Proportion of Cells") +
  ggtitle("Proportion of Cells in Liu et al.")

toplot <- bulk_eset@phenoData@data %>%
  dplyr::filter(group == "IRI_12m")
toplot <- toplot[5:ncol(toplot)]
toplot <- melt(toplot)
p4 <- ggplot(toplot, aes(x=variable, y=value, fill=variable)) + 
  geom_boxplot() + 
  ylab("Proportion of Cells") +
  ggtitle("Proportion of Cells in Liu et al.")

CombinePlots(list(p3,p4))

pdf("liu_deconvolution.pdf")
CombinePlots(list(p1,p2), nrow=2)
p3
dev.off()

