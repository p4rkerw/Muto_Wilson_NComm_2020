# this script will deconvolute the bulk RNAseq data for human IRI AKI in transplant patients
# data is obtained from PMID30429361
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
library(Biobase)

# read in the bulk eset using either
bulk_eset = readRDS("comparison_human_datasets/cippa_PMID30429361/eset.cippa.bulkRna.rds")

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

# calculate mean of PT_VCAM1 proportion across time points
pre <- pData(bulk_eset)[pData(bulk_eset)$group == "pre",]
mean(pre$PT_VCAM1)

post <- pData(bulk_eset)[pData(bulk_eset)$group == "post",]
mean(post$PT_VCAM1)

three <- pData(bulk_eset)[pData(bulk_eset)$group == "3months",]
mean(three$PT_VCAM1)

year <- pData(bulk_eset)[pData(bulk_eset)$group == "1year",]
mean(year$PT_VCAM1)

# determine if there is a statistically significant difference between mean PT_VCAM1
# proportion across time points using one-way ANOVA
toplot <- bulk_eset@phenoData@data
res.aov <- aov(PT_VCAM1 ~ group, data=toplot)
summary(res.aov) # p=9.7e-05

# compare pre- to post at 1 year
res.t <- t.test(pre$PT_VCAM1, year$PT_VCAM1)

# metadata located here
toplot <- bulk_eset@phenoData@data %>%
  dplyr::select(group, PT)
p1 <- ggplot(toplot, aes(x=group, y=PT, fill=group)) + 
  geom_boxplot() +
  ggtitle("Estimated Proportion of PT \nKidney Allograft IRI (Cippa 2018) \n(n=42)") +
  xlab("Time Point Post-Perfusion")

# plot proportion of cells in control group
toplot <- bulk_eset@phenoData@data %>%
  dplyr::select(group, PT_VCAM1)
p2 <- ggplot(toplot, aes(x=group, y=PT_VCAM1, fill=group)) + 
  geom_boxplot() +
  ggtitle("Estimated Proportion of PT_VCAM1 \nKidney Allograft IRI (Cippa 2018) \n(n=42)") +
  xlab("Time Point Post-Perfusion")


CombinePlots(list(p1,p2))

toplot <- bulk_eset@phenoData@data
toplot <- toplot[7:ncol(toplot)]
toplot <- melt(toplot)
ggplot(toplot, aes(x=variable, y=value, fill=variable)) + 
  geom_boxplot() + 
  ylab("Proportion of Cells") +
  ggtitle("Proportion of Cells")

pdf("cippa_deconvolution.pdf")
CombinePlots(list(p1,p2))
p3
dev.off()

