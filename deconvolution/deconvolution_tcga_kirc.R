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
sc.markers <- FindAllMarkers(rnaAggr) # only.pos=TRUE
list.markers <- sapply(levels(rnaAggr), function(x) {
  df <- dplyr::filter(sc.markers, cluster == x)  %>%
    # dplyr::slice(1:500)  %>%
    dplyr::select(gene)
})
names(list.markers) <- levels(rnaAggr)

# build reference basis matrix of expression profiles for individual celltypes
# average counts computed within each cell type in each sample
plotCellTotals(scrna_eset, 'cellType', 'SubjectName')
B <- bseqsc_basis(scrna_eset, list.markers, clusters = 'celltype', samples = 'SubjectName', ct.scale = TRUE)
plotBasis(B, list.markers, Colv = NA, Rowv = NA, layout = '_', col = 'Blues')

# read in bulk eset
bulk_eset <- readRDS("comparison_human_datasets/tcga_kirc/eset.tcga.bulkRna.rds")

# estimate celltype proportions with cibersort
fit <- bseqsc_proportions(bulk_eset, B, verbose = TRUE)

# add estimated celltype proportions as phenoData
pData(bulk_eset) <- cbind(pData(bulk_eset), t(coef(fit)))

# metadata located here
toplot <- bulk_eset@phenoData@data %>%
  dplyr::select(group, PT_KIM1)
p1 <- ggplot(toplot, aes(x=group, y=PT_KIM1, fill=group)) + 
  geom_boxplot() + 
  xlab("Patient Age Range") +
  ylab("Proportion of PT_VCAM1 Cells") +
  ggtitle("Proportion of PT_VCAM1 Cells by Patient Age Range TCGA-KIRC")

toplot$TCGA <- "TCGA"
p2 <- ggplot(toplot, aes(x=TCGA, y=PT_KIM1, fill=TCGA)) + 
  geom_boxplot() + 
  xlab("Patient Age Range") +
  ylab("Proportion of PT_VCAM1 Cells") +
  ggtitle("Proportion of PT_VCAM1 Cells by Patient Age Range TCGA-KIRC")


