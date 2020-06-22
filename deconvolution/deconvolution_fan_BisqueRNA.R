# devtools::install_github('shenorrlab/bseqsc')
library(here)
library(BisqueRNA) # for Seurat to ExpressionSet
library(Seurat)
library(dplyr)
library(ggplot2)
library(reshape)
library(forcats)

# read in bulk eset for diabetic nephropathy from fan et al
bulk_eset <- readRDS("comparison_human_datasets/fan_PMID31578193/eset.bulkRna.rds")

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
pData(bulk_eset)$group <- dplyr::recode(pData(bulk_eset)$group, A="Advanced", N="Control", B="Early")
pData(bulk_eset) <- dplyr::mutate(pData(bulk_eset), group = factor(group, levels = c("Control","Early","Advanced")))

# calculate mean of PT_VCAM1 proportion across groups
ctrl <- pData(bulk_eset)[pData(bulk_eset)$group == "Control",]
mean(ctrl$PT_VCAM1)

early <- pData(bulk_eset)[pData(bulk_eset)$group == "Early",]
mean(early$PT_VCAM1)

adv <- pData(bulk_eset)[pData(bulk_eset)$group == "Advanced",]
mean(adv$PT_VCAM1)
sd(adv$PT_VCAM1)

# determine if there is a statistically significant difference between mean PT_VCAM1
# proportion across groups using one-way ANOVA
toplot <- bulk_eset@phenoData@data
res.aov <- aov(PT_VCAM1 ~ group, data=toplot)
summary(res.aov) # p=0.71

# metadata located here
# group A=advanced, B=early, N=normal 
toplot <- bulk_eset@phenoData@data %>%
  dplyr::select(group, PT) 
p1 <- ggplot(toplot, aes(x=group, y=PT, fill=group)) + 
  xlab("") +
  geom_boxplot() +
  ggtitle("Estimated Proportion of PT \nHuman Diabetic Nephropathy \n(Fan 2019) (n=36)") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

toplot <- bulk_eset@phenoData@data %>%
  dplyr::select(group, PT_VCAM1) 
p2 <- ggplot(toplot, aes(x=group, y=PT_VCAM1, fill=group)) + 
  xlab("") +
  geom_boxplot() +
  ggtitle("Estimated Proportion of PT_VCAM1 \nHuman Diabetic Nephropathy \n(Fan 2019) (n=36)")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


CombinePlots(list(p1,p2))

# plot proportion of all cell types across groups
toplot <- bulk_eset@phenoData@data
toplot <- dplyr::select(toplot, -lib.size, -norm.factors, -Run, -GEO_Accession,
                        -Sample_ID, -replaceable)
toplot <- melt(toplot)
p3 <- ggplot(toplot, aes(x=variable, y=value, fill=variable)) + 
  geom_boxplot() + 
  ylab("Proportion of Cells") +
  ggtitle("Proportion of Cells in all Groups") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


pdf("diabetes_deconvolution.pdf")
CombinePlots(list(p1,p2), nrow=2)
dev.off()





# adjust the degs pval for celltype proportions
fit_edger <- fitEdgeR(bulk_eset, ~group, coef = c("groupA","groupB"))

write.csv(fit_edger, file = "comparison_human_datasets/fan_PMID31578193/adjusted_deg.csv")






