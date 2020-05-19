# devtools::install_github('shenorrlab/bseqsc')
library(bseqsc)
library(here)
library(BisqueRNA) # for Seurat to ExpressionSet
library(Seurat)
library(dplyr)
library(ggplot2)
library(reshape)
bseqsc_config('comparison_human_datasets/fan_PMID31578193/CIBERSORT.R')


rnaAggr <- readRDS("cellranger_rna_prep/rnaAggr_control.rds")
scrna_eset <- SeuratToExpressionSet(rnaAggr, delimiter="-", position=2, version="v3")
scrna_eset@phenoData$celltype <- rnaAggr$celltype

# select top cluster specific snRNA markers for deconvolution. 
sc.markers <- FindAllMarkers(rnaAggr, min.pct = 0.25, only.pos = TRUE) # only.pos=TRUE
list.markers <- sapply(levels(rnaAggr), function(x) {
  df <- dplyr::filter(sc.markers, cluster == x)  %>%
    dplyr::filter(avg_logFC > 1 & p_val_adj < 0.01) %>%
    # dplyr::arrange(desc(avg_logFC, p_val_adj)) %>%
    # dplyr::slice(1:150)  %>%
    dplyr::select(gene)
})
names(list.markers) <- levels(rnaAggr)

# build reference basis matrix of expression profiles for individual celltypes
# average counts computed within each cell type in each sample
plotCellTotals(scrna_eset, 'cellType', 'SubjectName')
B <- bseqsc_basis(scrna_eset, list.markers, clusters = 'celltype', samples = 'SubjectName', ct.scale = TRUE)
plotBasis(B, list.markers, Colv = NA, Rowv = NA, layout = '_', col = 'Blues')

# load the bulk eset
bulk_eset <- readRDS("comparison_mouse_datasets/liu_PMID28931758/eset.liu.bulkRna.rds")

# estimate celltype proportions with cibersort
fit <- bseqsc_proportions(bulk_eset, B, verbose = TRUE)

# add estimated celltype proportions as phenoData
pData(bulk_eset) <- cbind(pData(bulk_eset), t(coef(fit)))

# metadata located here
bulk_eset@phenoData@data$group <- factor(bulk_eset@phenoData@data$group, 
                                            levels=c("NORM3m","NORM9m","NORM15m","SHAM12m","SHAM4h","SHAM24h",
                                                     "IRI6mN","IRI2h","IRI4h","IRI24h","IRI48h","IRI72h",
                                                     "IRI7d","IRI14d","IRI28d","IRI12m"))

toplot <- bulk_eset@phenoData@data %>%
  dplyr::select(group, PT_VCAM1)
p1 <- ggplot(toplot, aes(x=group, y=PT_VCAM1, fill=group)) + 
  geom_boxplot() + 
  xlab("Treatment Group") +
  ylab("Proportion of PT_VCAM1 Cells") +
  ggtitle("Proportion of PT_VCAM1 Cells in Murine IRI") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  NoLegend()


toplot <- bulk_eset@phenoData@data %>%
  dplyr::select(group, PT)
p2 <- ggplot(toplot, aes(x=group, y=PT, fill=group)) + 
  geom_boxplot() + 
  xlab("Treatment Group") +
  ylab("Proportion of PT Cells") +
  ggtitle("Proportion of PT Cells in Murine IRI") +
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

toplot <- bulk_eset@phenoData@data
toplot <- toplot[5:ncol(toplot)]
toplot <- melt(toplot)
p3 <- ggplot(toplot, aes(x=variable, y=value, fill=variable)) + 
  geom_boxplot() + 
  ylab("Proportion of Cells") +
  ggtitle("Proportion of Cells in Liu et al.")

toplot <- bulk_eset@phenoData@data
toplot$aggregated <- str_extract(toplot$group, pattern = "[A-Z]+")
toplot <-
  dplyr::select(toplot, aggregated, PT_KIM1)
p4 <- ggplot(toplot, aes(x=aggregated, y=PT_KIM1, fill=aggregated)) + 
  geom_boxplot() + 
  xlab("Treatment Group") +
  ylab("Proportion of PT_VCAM1 Cells") +
  ggtitle("Proportion of PT_VCAM1 Cells in Murine IRI") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  NoLegend()


pdf("liu_deconvolution.pdf")
CombinePlots(list(p1,p2), nrow=2)
p3
dev.off()

