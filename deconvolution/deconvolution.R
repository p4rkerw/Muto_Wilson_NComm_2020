# this script will deconvolute the bulk RNAseq data for human IRI AKI in transplant patients
# data is obtained from PMID30429361

# devtools::install_github('shenorrlab/bseqsc')
library(bseqsc)
library(here)
library(BisqueRNA) # for Seurat to ExpressionSet
library(Seurat)
library(dplyr)
library(ggplot2)
library(reshape)
library(data.table)
library(DEFormats)
library(tibble)
library(tximport)
library(biomaRt)
bseqsc_config('github_repository/healthyKidney/deconvolution/CIBERSORT.R')

bulk_eset = readRDS("comparison_human_datasets/cippa_PMID30429361/eset.cippa.bulkRna.rds")

# load the snrna data
rnaAggr <- readRDS("cellranger_rna_prep/rnaAggr_control.rds")

# convert the gene names to ensembl id for the snrna data
scrna_eset <- SeuratToExpressionSet(rnaAggr, delimiter="-", position=2, version="v3")
scrna_eset@phenoData$celltype <- rnaAggr$celltype

# select top 50 cluster specific snRNA markers for deconvolution. 
# downsample to 1000 cells per cluster to speed up computation
sc.markers <- FindAllMarkers(rnaAggr, min.pct = 0.25, max.cells.per.ident=1000, only.pos = TRUE)
list.markers <- sapply(levels(rnaAggr), function(x) {
  df <- dplyr::filter(sc.markers, cluster == x)  %>%
    dplyr::slice(1:100)  %>%
    dplyr::select(gene)
})
names(list.markers) <- levels(rnaAggr)

# build reference basis matrix of expression profiles for individual celltypes
# average counts computed within each cell type in each sample
plotCellTotals(scrna_eset, 'cellType', 'SubjectName')
B <- bseqsc_basis(scrna_eset, list.markers, clusters = 'celltype', 
                  samples = 'SubjectName', ct.scale = TRUE)
plotBasis(B, list.markers, Colv = NA, Rowv = NA, layout = '_', col = 'Blues')


# estimate celltype proportions with cibersort
fit <- bseqsc_proportions(bulk_eset, B, verbose = TRUE)

# add estimated celltype proportions as phenoData
pData(bulk_eset) <- cbind(pData(bulk_eset), t(coef(fit)))

# metadata located here
toplot <- bulk_eset@phenoData@data %>%
  dplyr::select(group, PT)
p1 <- ggplot(toplot, aes(x=group, y=PT, fill=group)) + 
  geom_boxplot() +
  ggtitle("Proportion of PT (Cippa 2018)")



# plot proportion of cells in control group
toplot <- bulk_eset@phenoData@data %>%
  dplyr::select(group, PT_KIM1)
p2 <- ggplot(toplot, aes(x=group, y=PT_KIM1, fill=group)) + 
  geom_boxplot() +
  ggtitle("Proportion of PT_VCAM1 (Cippa 2018)")


CombinePlots(list(p1,p2))

# adjust the degs pval for celltype proportions
fit_edger <- fitEdgeR(bulk_eset, ~group, coef = c("post","3months", "1year"))

#write.csv(fit_edger, file = "comparison_human_datasets/cippa_PMID30429361/adjusted_deg.csv")






