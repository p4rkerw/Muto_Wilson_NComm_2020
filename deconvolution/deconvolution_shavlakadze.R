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


# ALTERNATE: select cluster markers
# celltype.markers <- c("CUBN","HAVCR1","SLC5A1","SLC5A2", # PT and PT-KIM1+ markers
#                       "CFH", # PEC
#                       "SLC12A1", # TAL NKCC2
#                       "SLC12A3","TRPM6", # DCT1 and DCT2 NCC
#                       "SCNN1G","TRPV5", # DCT2/CNT ENaC
#                       "CALB1", # CNT
#                       "AQP2", # PC
#                       "ATP6V0D2", # ICA and ICB
#                       "SLC4A1","SLC26A7", # ICA
#                       "SLC26A4", # ICB
#                       "NPHS1","NPHS2", # PODO
#                       "PECAM1","FLT1","EMCN", # ENDO
#                       "CLDN5", # GEC
#                       "ITGA8","PDGFRB", # MES
#                       "ACTA2","CALD1", # FIB
#                       "PTPRC") # WBC


# build reference basis matrix of expression profiles for individual celltypes
# average counts computed within each cell type in each sample
plotCellTotals(scrna_eset, 'cellType', 'SubjectName')
B <- bseqsc_basis(scrna_eset, list.markers, clusters = 'celltype', samples = 'SubjectName', ct.scale = TRUE)
plotBasis(B, list.markers, Colv = NA, Rowv = NA, layout = '_', col = 'Blues')

# read in bulk eset
bulk_eset <- readRDS("comparison_mouse_datasets/shavlakadze_PMID31533046/eset.shavlakadze.bulkRna.rds")

# estimate celltype proportions with cibersort
fit <- bseqsc_proportions(bulk_eset, B, verbose = TRUE)

# add estimated celltype proportions as phenoData
pData(bulk_eset) <- cbind(pData(bulk_eset), t(coef(fit)))

# metadata located here
bulk_eset@phenoData@data$group <- factor(bulk_eset@phenoData@data$group, levels=c("6m","9m","12m","18m","21m","24m","27m"))

toplot <- bulk_eset@phenoData@data %>%
  dplyr::select(group, PT_KIM1)
p1 <- ggplot(toplot, aes(x=group, y=PT_KIM1, fill=group)) + 
  geom_boxplot() + 
  xlab("Rat Age") +
  ylab("Proportion of PT_VCAM1 Cells") +
  ggtitle("Proportion of PT_VCAM1 Cells by Rat Age")

toplot <- bulk_eset@phenoData@data %>%
  dplyr::select(group, PT)
p2 <- ggplot(toplot, aes(x=group, y=PT, fill=group)) + 
  geom_boxplot() + 
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


pdf("shavlakadze_deconvolution.pdf")
CombinePlots(list(p1,p2))
p3
dev.off()