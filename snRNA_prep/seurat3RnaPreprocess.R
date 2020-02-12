# this script will preprocess, filter, and annotate aggregated healthy kidney snRNAseq libraries
library(Seurat) # 3.0.2
library(ggplot2)
library(harmony) # 1.0
library(here) # set all filepaths relative to analysis directory of an R project eg. G:/diabneph
set.seed(1234)
sessionInfo()
here()

outs <- "cellranger_rna_aggr_control/outs/"

# load aggregated snRNAseq data from cellranger aggregate matrix and create seurat objects
counts <- Read10X(here(outs,"filtered_feature_bc_matrix"))
metadata <- read.csv(here(outs,"aggregation.csv"))

# Not implemented: filter out barcodes identified by doubletfinder as doublets

# create seurat object
rnaAggr <- CreateSeuratObject(counts = counts, min.cells = 10, min.features = 500, project = "RNA")

# extract GEM groups from individual barcodes using string split and the suffix integer
# use the GEM groups to assign sample origin (control vs. diabetes) from the aggregation.csv metadata
# files were aggregated control 1, 2, 3, diabetes 1, 2, 3 are correspond to GEM groups 1-6
gemgroup <- sapply(strsplit(rownames(rnaAggr@meta.data), split="-"), "[[", 2) 
current.gemgroups <- seq(1, length(levels(metadata$library_id)))
orig.ident <- levels(metadata$library_id) 
sampleID <- plyr::mapvalues(gemgroup, from = current.gemgroups, to = orig.ident)
rnaAggr <- AddMetaData(object=rnaAggr, metadata=data.frame(orig.ident=sampleID, row.names=rownames(rnaAggr@meta.data)))
rnaAggr <- PercentageFeatureSet(rnaAggr, pattern = "^MT-", col.name = "percent.mt")
rnaAggr <- PercentageFeatureSet(rnaAggr, pattern = "^RPL", col.name = "percent.rpl")
rnaAggr <- PercentageFeatureSet(rnaAggr, pattern = "^RPS", col.name = "percent.rps")
# control vs diabetes
current.ids <- orig.ident
new.ids <- metadata$group
rnaAggr@meta.data$disease<-plyr::mapvalues(rnaAggr@meta.data$orig.ident,from = current.ids,to=new.ids)

# visualize filtering parameters
# VlnPlot(rnaAggr, features = c("nCount_RNA","nFeature_RNA","percent.mt"), pt.size = 0.1)

# Filtering low quality cells with high mitochondrial RNAs and ribosomal protein RNAs (known to be stable and thus enriched in cells compared to nuclei)
rnaAggr <- subset(rnaAggr, 
                  subset = nFeature_RNA > 500 
                  & nFeature_RNA < 4000 
                  & nCount_RNA < 16000 
                  & percent.mt < 0.8 
                  & percent.rps < 0.4 
                  & percent.rpl < 0.4)

# run sctransform
# Regress out the mitochondrial reads and nCount_RNA
rnaAggr <- SCTransform(rnaAggr, vars.to.regress = c("percent.mt","percent.rpl", "percent.rps","nCount_RNA"), verbose = TRUE)
rnaAggr <- RunPCA(rnaAggr, verbose = TRUE)
# ElbowPlot(rnaAggr, ndims = 50) # to determin number of dimensions for clustering
rnaAggr <- RunHarmony(rnaAggr, "orig.ident", plot_convergence = TRUE)
rnaAggr <- FindNeighbors(rnaAggr, dims = 1:30, verbose = TRUE, reduction = "harmony")
rnaAggr <- FindClusters(rnaAggr, verbose = TRUE, resolution = 0.6, reduction = "harmony")
rnaAggr <- RunUMAP(rnaAggr, dims = 1:30, verbose = TRUE, reduction = "harmony")

# visualize the clustering
# DimPlot(rnaAggr, reduction = "UMAP", assay = "SCT")
p1 <- DimPlot(rnaAggr, reduction = "umap", assay = "SCT", label = TRUE) + ggtitle("snRNA Original Clustering")

celltype.markers <- c("CUBN","LRP2","HAVCR1","SLC5A1","SLC5A2", # PT and PT-KIM1+ markers
                      "CFH", # PEC
                      "SLC12A1", # TAL NKCC2
                      "SLC12A3","TRPM6", # DCT1 and DCT2 NCC
                      "SCNN1G","TRPV5", # DCT2/CNT ENaC
                      "CALB1", # CNT
                      "AQP2", # PC
                      "ATP6V0D2", # ICA and ICB
                      "SLC4A1", # ICA
                      "SLC26A4", # ICB
                      "NPHS1","NPHS2", # PODO
                      "PECAM1","FLT1", # ENDO
                      "CLDN5", # GEC
                      "ITGA8","PDGFRB", # MES
                      "PTPRC") # WBC
p2 <- DotPlot(rnaAggr, features = celltype.markers) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

CombinePlots(plots = list(p1, p2))
########## Annotation of the clusters for data integration with snATAC dataset ##########
rnaAggr[["orig.clusters"]] <- Idents(object = rnaAggr) # stash cluster idents prior to annotation
new.cluster.ids <- c("PCT","TAL","DCT1","CNT","PCT","ICA","TAL","PC","PEC","ENDO","DCT2","PT_KIM1","PODO","MES","ICB","TAL","ENDO","ENDO","LEUK")
names(new.cluster.ids) <- levels(rnaAggr)
rnaAggr <- RenameIdents(rnaAggr, new.cluster.ids)

# reorder the idents and save celltype annotations in celltype slot metadata
levels(rnaAggr) <- c("PCT","PT_KIM1","PEC","TAL","DCT1","DCT2","CNT","PC","ICA","ICB","PODO","ENDO","MES","LEUK")
rnaAggr@meta.data$celltype <- rnaAggr@active.ident

# create low-resolution celltype identities for snATAC thresholding (ie group PT and PT-KIM1 and distal nephron together)
lowres.cluster.ids <- c("PCT","PCT","PEC","TAL","DCT_CNT_PC","DCT_CNT_PC","DCT_CNT_PC","DCT_CNT_PC",
                        "ICA","ICB","PODO","ENDO","MES","LEUK")
names(lowres.cluster.ids) <- levels(rnaAggr)
rnaAggr <- RenameIdents(rnaAggr, lowres.cluster.ids)
rnaAggr[["lowres.celltype"]] <- Idents(rnaAggr)

# reset celltype as primary ident before saving and plotting
Idents(rnaAggr) <- "celltype"

# redraw umap and dotplot with reordered idents
p3 <- DimPlot(rnaAggr, reduction = "umap", assay = "SCT", label = TRUE) + ggtitle("snRNA Annotated Celltypes")
p4 <- DotPlot(rnaAggr, features = celltype.markers) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

# draw pdf plots for before and after annotation
dir.create("plots", showWarnings = FALSE)
pdf(here("plots","umap.rnaAggr.pdf"))
p1
p2
p3
p4
dev.off()

print("Saving aggregated snRNAseq object as rnaAggr.rds in:")
dir.create("cellranger_rna_prep", showWarnings = FALSE)
here("cellranger_rna_prep")
saveRDS(rnaAggr, file = here("cellranger_rna_prep","rnaAggr_control.rds"))
        