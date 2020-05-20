# devtools::install_github('shenorrlab/bseqsc')
library(bseqsc)
library(here)
library(BisqueRNA) # for Seurat to ExpressionSet
library(Seurat)
library(dplyr)
library(ggplot2)
library(reshape)
bseqsc_config('comparison_human_datasets/fan_PMID31578193/CIBERSORT.R')

bulk_eset <- readRDS("comparison_human_datasets/fan_PMID31578193/eset.bulkRna.rds")

rnaAggr <- readRDS("cellranger_rna_prep/rnaAggr_combined.rds")
scrna_eset <- SeuratToExpressionSet(rnaAggr, delimiter="-", position=2, version="v3")
scrna_eset@phenoData$celltype <- rnaAggr$celltype

# select top 50 cluster specific snRNA markers for deconvolution. 
# downsample to 1000 cells per cluster to speed up computation
sc.markers <- FindAllMarkers(rnaAggr, min.pct = 0.25, max.cells.per.ident=1000, only.pos = TRUE)
list.markers <- sapply(levels(rnaAggr), function(x) {
  df <- dplyr::filter(sc.markers, cluster == x)  %>%
    dplyr::slice(1:50)  %>%
    dplyr::select(gene)
})
names(list.markers) <- levels(rnaAggr)

# build reference basis matrix of expression profiles for individual celltypes
# average counts computed within each cell type in each sample
plotCellTotals(scrna_eset, 'cellType', 'SubjectName')
B <- bseqsc_basis(scrna_eset, list.markers, clusters = 'celltype', samples = 'SubjectName', ct.scale = TRUE)
plotBasis(B, list.markers, Colv = NA, Rowv = NA, layout = '_', col = 'Blues')


# estimate celltype proportions with cibersort
fit <- bseqsc_proportions(bulk_eset, B, verbose = TRUE)

# add estimated celltype proportions as phenoData
pData(bulk_eset) <- cbind(pData(bulk_eset), t(coef(fit)))

# metadata located here
# group A=advanced, B=early, N=normal 
toplot <- bulk_eset@phenoData@data %>%
  dplyr::select(group, PT) %>%
  dplyr::mutate(group = factor(group, labels=c("Control","Advanced","Early")))
ggplot(toplot, aes(x=group, y=PT, fill=group)) + 
  geom_boxplot()


# plot proportion of cells in control group
toplot <- bulk_eset@phenoData@data %>%
  dplyr::mutate(group = factor(group, labels=c("Control","Advanced","Early"))) %>%
  dplyr::filter(group == "Control") %>%
  dplyr::select(-replaceable)
toplot <- toplot[, 6:ncol(toplot)] # select group id and celltype proportions
toplot2 <- melt(toplot)
ggplot(toplot2, aes(x=variable, y=value)) +
  ggtitle("Celltype proportion in bulk RNA-seq of Healthy Control Kidney")
  geom_boxplot()

# plot proportion of cells in advanced group
toplot <- bulk_eset@phenoData@data %>%
  dplyr::mutate(group = factor(group, labels=c("Control","Advanced","Early"))) %>%
  dplyr::filter(group == "Advanced") %>%
  dplyr::select(-replaceable)
toplot <- toplot[, 6:ncol(toplot)] # select group id and celltype proportions
toplot2 <- melt(toplot)
ggplot(toplot2, aes(x=variable, y=value)) + 
  ggtitle("Celltype proportion in bulk RNA-seq of Advanced Diabetic Kidney") +
  geom_boxplot()

# adjust the degs pval for celltype proportions
fit_edger <- fitEdgeR(bulk_eset, ~group, coef = c("groupA","groupB"))

write.csv(fit_edger, file = "comparison_human_datasets/fan_PMID31578193/adjusted_deg.csv")






