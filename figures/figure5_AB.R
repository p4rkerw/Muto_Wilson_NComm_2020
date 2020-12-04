# load the gene count files from Liu et al PMID28931758 
library(readr)
library(dplyr)
library(stringr)
library(data.table)
library(DESeq2)
library(biomaRt)
library(tibble)
library(GEOquery)
library(purrr)
library(openxlsx)
# library(EnsDb.Hsapiens.v86)

# read in the FPKM counts 
mouseCounts <- read.xlsx("comparison_mouse_datasets/liu_PMID28931758/GSE98622_mouse-iri-master.xlsx") %>%
  dplyr::filter(grepl("ENSM", X1)) %>%
  dplyr::distinct(symbol, .keep_all = TRUE)

# convert mouse genes to human genes
mouseGenes <- mouseCounts$symbol
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
lookup_table = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = mouseGenes , mart = mouse,
                      attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)

# extract the counts from genes that are shared between species, arrange by mean count per transcript and keep
# the transcript with the most counts
colnames(lookup_table)[1] <- "symbol"
counts <- merge(x=lookup_table, y=mouseCounts, by="symbol") 
counts$mean <- rowMeans(as.matrix(counts[5:ncol(counts)]), na.rm = TRUE)
counts <- dplyr::group_by(counts, HGNC.symbol) %>%
  dplyr::arrange(desc(mean), .by_group = TRUE) %>%
  dplyr::distinct(HGNC.symbol, .keep_all=TRUE) %>%
  dplyr::select(-X1, -type, -symbol, -mean) %>%
  column_to_rownames(var="HGNC.symbol")

# convert to integer counts
i <- seq(colnames(counts))
counts[ , i] <- apply(counts[ , i], 2,            
                      function(x) as.integer(x))

# generate coldata
coldata <- data.frame(sample = colnames(counts))
coldata$condition <- str_split(coldata$sample, pattern = "-", simplify = TRUE)[,1]

# import data into deseq2 for rna seq analysis and level the conditions
dds <- DESeqDataSetFromMatrix(counts,
                              colData = coldata,
                              design = ~ condition)

# do not put through DESeq2 because these are FPKM as opposed to read counts or transcript abundance
# dds <- DESeq(dds)
# counts <- counts(dds, normalized=TRUE)[-1]
# do not normalize the FPKM values

# convert dds to a dgelist object and then an expressionset to prepare for deconvolution 
library(DEFormats)
dge <- as.DGEList(dds)
pheno <- as(dge$samples, "AnnotatedDataFrame") # coerces a data.frame to an AnnotatedDataFrame.
eset <- ExpressionSet(assayData = dge$counts, phenoData = pheno)
saveRDS(eset, file = "comparison_mouse_datasets/liu_PMID28931758/eset.liu.bulkRna.rds")

###################################################################################################################

# devtools::install_github('shenorrlab/bseqsc')
library(bseqsc)
library(here)
library(BisqueRNA) # for Seurat to ExpressionSet
library(Seurat)
library(dplyr)
library(ggplot2)
library(reshape)
library(ggpubr)
library(multcomp)
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
bulk_eset_all <- readRDS("comparison_mouse_datasets/liu_PMID28931758/eset.liu.bulkRna.rds")

# add estimated celltype proportions as phenoData
res <- BisqueRNA::ReferenceBasedDecomposition(bulk_eset_all, scrna_eset, markers=NULL, use.overlap=F)
pData(bulk_eset_all) <- cbind(pData(bulk_eset_all), t(res[["bulk.props"]]))

# metadata located here
bulk_eset_all@phenoData@data$group <- factor(bulk_eset_all@phenoData@data$group, 
                                         levels=c("SHAM4h","SHAM24h","SHAM12m",
                                                  "IRI2h","IRI4h","IRI24h","IRI48h","IRI72h",
                                                  "IRI7d","IRI14d","IRI28d","IRI6m","IRI12m","NORM3m","NORM9m","NORM15m"))
#rownames(bulk_eset_all@phenoData)
#[1] "SHAM4h-R1"  "SHAM4h-R2"  "SHAM4h-R3"  "SHAM24h-R1" "SHAM24h-R2" "SHAM24h-R3" "NORM3m-R1" 
#[8] "NORM3m-R2"  "NORM3m-R3"  "NORM9m-R1"  "NORM9m-R2"  "NORM9m-R3"  "SHAM12m-R1" "SHAM12m-R2"
#[15] "SHAM12m-R3" "NORM15m-R1" "NORM15m-R2" "NORM15m-R3" "IRI2h-R1"   "IRI2h-R2"   "IRI2h-R3"  
#[22] "IRI4h-R1"   "IRI4h-R2"   "IRI4h-R3"   "IRI24h-R1"  "IRI24h-R2"  "IRI24h-R3"  "IRI48h-R1" 
#[29] "IRI48h-R2"  "IRI48h-R3"  "IRI72h-R1"  "IRI72h-R2"  "IRI72h-R3"  "IRI7d-R1"   "IRI7d-R2"  
#[36] "IRI7d-R3"   "IRI14d-R1"  "IRI14d-R2"  "IRI14d-R3"  "IRI28d-R1"  "IRI28d-R2"  "IRI28d-R3" 
#[43] "IRI6m-R1"   "IRI6m-R2"   "IRI6m-R3"   "IRI6m-R4"   "IRI12m-R1"  "IRI12m-R2"  "IRI12m-R3" 


###########IRI############### 

bulk_eset <- bulk_eset_all[,c(1:6,13:15,19:49)] #IRI dataset

toplot <- bulk_eset@phenoData@data %>%
  dplyr::select(group, c(PT,PT_VCAM1,LEUK))

fig5a_1 <- ggboxplot(toplot, x = "group", y = "PT_VCAM1",
          add = "jitter",ylim = c(0, 0.075),
          fill = "group", palette = c("#CCDFE7","#CCDFE7","#CCDFE7","#FFFDDE",jdb_palette("brewer_fire")))+NoLegend()  #590x480

fig5a_2 <- ggboxplot(toplot, x = "group", y = "PT",
                   add = "jitter",ylim = c(0, 0.4),
                   fill = "group", palette = c("#CCDFE7","#CCDFE7","#CCDFE7","#FFFDDE",jdb_palette("brewer_fire")))+NoLegend() 



#p value
bulk_eset_all@phenoData@data$group <- factor(bulk_eset_all@phenoData@data$group, 
                                          levels=c("SHAM24h","SHAM4h","SHAM12m",
                                                   "IRI2h","IRI4h","IRI24h","IRI48h","IRI72h",
                                                   "IRI7d","IRI14d","IRI28d","IRI6m","IRI12m","NORM3m","NORM9m","NORM15m"))
#control is IRI24h
bulk_eset <- bulk_eset_all[,c(1:6,13:15,19:49)] 
toplot <- bulk_eset@phenoData@data %>%
  dplyr::select(group, c(PT,PT_VCAM1,LEUK))
res.aov <- aov(PT_VCAM1 ~ group, data =toplot)
summary(res.aov)
# Summary of the analysis
glht.res <- glht(res.aov, linfct = mcp(group = "Dunnett"))
summary(glht.res)

###########Aging############### 

bulk_eset <- bulk_eset_all[,c(7:12,16:18)] #aging  

toplot <- bulk_eset@phenoData@data %>%
  dplyr::select(group, c(PT,PT_VCAM1,LEUK))

#290x480

fig5b_1 <- ggboxplot(toplot, x = "group", y = "PT_VCAM1",
                     add = "jitter",ylim = c(0, 0.075),
                     fill = "group", palette = c("#DABC77","#AB8939","#674E14",jdb_palette("brewer_fire")))+NoLegend() 

fig5b_2 <- ggboxplot(toplot, x = "group", y = "PT",
                     add = "jitter",ylim = c(0, 0.4),
                     fill = "group", palette = c("#DABC77","#AB8939","#674E14",jdb_palette("brewer_fire")))+NoLegend() 

fig5b_3 <- ggboxplot(toplot, x = "group", y = "LEUK",
                     add = "jitter",ylim = c(0, 0.075),
                     fill = "group", palette = c("#DABC77","#AB8939","#674E14",jdb_palette("brewer_fire")))+NoLegend() 

#p value
res.aov <- aov(PT_VCAM1 ~ group, data =toplot)
summary(res.aov)

# Summary of the analysis
glht.res <- glht(res.aov, linfct = mcp(group = "Dunnett"))
summary(glht.res)

