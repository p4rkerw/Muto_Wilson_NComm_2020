library(here)
library(BisqueRNA) # for Seurat to ExpressionSet
library(Seurat)
library(dplyr)
library(ggplot2)
library(reshape)
library(Biobase)

# read in bulk eset
bulk_eset <- readRDS("comparison_human_datasets/tcga_kirc/eset.tcga.bulkRna.rds")

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
age <- pData(bulk_eset)$age
age_group <- cut(age, breaks=c(30,60,70,80,90))
pData(bulk_eset)$age_group <- age_group

# metadata located here
toplot <- bulk_eset@phenoData@data %>%
  dplyr::select(age_group, age, PT_VCAM1)

# determine if there is a statistically significant difference between mean PT_VCAM1
# proportion across groups using one-way ANOVA
res.aov <- aov(PT_VCAM1 ~ age_group, data=toplot)
summary(res.aov) # p=0.71

# visualize a potential linear relationship between age and PT_VCAM1 proportion
scatter.smooth(x=toplot$age, y=toplot$PT_VCAM1, main="PT_VCAM1 vs Age")
cor(toplot$age, toplot$PT_VCAM1) # correlation is 0.14

# test for significant with a glm
linearMod <- lm(PT_VCAM1 ~ age, data=toplot)
summary(linearMod) # p=0.32

# metadata located here
toplot <- bulk_eset@phenoData@data %>%
  dplyr::select(group, PT_VCAM1)
p1 <- ggplot(toplot, aes(x=group, y=PT_VCAM1, fill=group)) + 
  geom_boxplot() + 
  xlab("Patient Age Range") +
  ylab("Proportion of PT_VCAM1 Cells") +
  ggtitle("Proportion of PT_VCAM1 Cells by Patient Age Range TCGA")

toplot <- bulk_eset@phenoData@data %>%
  dplyr::select(age_group, PT_VCAM1)
p2 <- ggplot(toplot, aes(x=age_group, y=PT_VCAM1, fill=age_group)) + 
  geom_boxplot() + 
  xlab("Patient Age Range") +
  ylab("Proportion of PT_VCAM1 Cells") +
  ggtitle("Proportion of PT_VCAM1 Cells by Patient Age Range TCGA")


toplot <- bulk_eset@phenoData@data
toplot <- dplyr::select(toplot, all_of(unique(new.cluster.ids)))
toplot <- melt(toplot)
p3 <- ggplot(toplot, aes(x=variable, y=value, fill=variable)) + 
  geom_boxplot() + 
  ylab("Proportion of Cells") +
  xlab("") +
  ggtitle("Estimated Proportion of Cells in Non-Tumor TCGA Bulk RNA-seq (n=72)") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


