# devtools::install_github('shenorrlab/bseqsc')
library(bseqsc)
library(here)
library(BisqueRNA) # for Seurat to ExpressionSet
library(Seurat)
library(dplyr)
library(ggplot2)
library(reshape)
library(ggpubr)
library(BuenColors)

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
bulk_eset <- readRDS("comparison_human_datasets/fan_PMID31578193/eset.fan.bulkRna.rds")

# add estimated celltype proportions as phenoData
res <- BisqueRNA::ReferenceBasedDecomposition(bulk_eset, scrna_eset, markers=NULL, use.overlap=F)
pData(bulk_eset) <- cbind(pData(bulk_eset), t(res[["bulk.props"]]))


# metadata located here
bulk_eset@phenoData@data$group <- factor(bulk_eset@phenoData@data$group, 
                                         levels=c("N","B","A"))

toplot <- bulk_eset@phenoData@data %>%
  dplyr::select(group, c(PT,PT_VCAM1,LEUK))
#330x480
fig5d_1 <- ggboxplot(toplot, x = "group", y = "PT_VCAM1",
          add = "jitter",ylim = c(0, 0.1),
          fill = "group", palette = c("#A3CDE5","#E8BD72","#D9A1F3"))+NoLegend() 

fig5d_2 <- ggboxplot(toplot, x = "group", y = "PT",
                   add = "jitter",ylim = c(0, 0.4),
                   fill = "group", palette = c("#A3CDE5","#E8BD72","#D9A1F3"))+NoLegend() 

fig5d_3 <- ggboxplot(toplot, x = "group", y = "LEUK",
                   add = "jitter",ylim = c(0, 0.075),
                   fill = "group", palette = c("#A3CDE5","#E8BD72","#D9A1F3"))+NoLegend() 

#toplot <- bulk_eset@phenoData@data
#toplot <- toplot[8:ncol(toplot)]
#toplot <- melt(toplot)
#fig5a_4 <- ggboxplot(toplot, x = "variable", y = "value",
                     #add = "jitter",ylim = c(0, 0.5),
                     #fill = "variable", palette = c("#CCDFE7","#C2318E","#CCDFE7",jdb_palette("brewer_marine")))+NoLegend() 
#882x480

#p value
library(multcomp)
res.aov <- aov(PT_VCAM1 ~ group, data =toplot)
summary(res.aov)
# Summary of the analysis
glht.res <- glht(res.aov, linfct = mcp(group = "Dunnett"))
summary(glht.res)