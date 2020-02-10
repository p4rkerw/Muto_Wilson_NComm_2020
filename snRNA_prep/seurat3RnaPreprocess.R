# this script will preprocess, filter, and annotate aggregated healthy kidney snRNAseq libraries
library(Seurat) # 3.0.2
library(ggplot2)
library(harmony) # 1.0
library(here) # set all filepaths relative to analysis directory eg. G:/diabneph
set.seed(1234)
sessionInfo()
here()

outs <- "cellranger_rna_counts/version_3.1.0/cellranger_rna_aggr_control/outs/"

# load aggregated snRNAseq data from cellranger aggregate matrix and create seurat objects
counts <- Read10X(here(outs,"filtered_feature_bc_matrix"))
metadata <- read.csv(here(outs,"aggregation.csv"))
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

# Filtering low quality cells with high mitochondrial RNAs and ribosomal protein RNAs (known to be stable and thus enriched in cells compared to nuclei)
rnaAggr <- subset(rnaAggr, 
                  subset = nFeature_RNA > 500 
                  & nFeature_RNA < 6000 
                  & nCount_RNA < 16000 
                  & percent.mt < 0.8 
                  & percent.rps < 0.4 
                  & percent.rpl < 0.4)

# run sctransform
# Regress out the mitochondrial reads and nCount_RNA
rnaAggr <- SCTransform(rnaAggr, vars.to.regress = c("percent.mt","percent.rpl", "percent.rps","nCount_RNA"), verbose = TRUE)
rnaAggr <- RunPCA(rnaAggr, verbose = TRUE)
# ElbowPlot(rnaAggr, ndims = 40) # to determin number of dimensions for clustering
rnaAggr <- RunHarmony(rnaAggr, "orig.ident", plot_convergence = TRUE)
rnaAggr <- FindNeighbors(rnaAggr, dims = 1:20, verbose = TRUE, reduction = "harmony")
rnaAggr <- FindClusters(rnaAggr, verbose = TRUE, resolution = 0.6, reduction = "harmony")
rnaAggr <- RunUMAP(rnaAggr, dims = 1:20, verbose = TRUE, reduction = "harmony")

# visualize the clustering
# DimPlot(rnaAggr, reduction = "UMAP", assay = "SCT")
p1 <- DimPlot(rnaAggr, reduction = "umap", assay = "SCT", label = TRUE)

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
                      "ITGA8","PDGFRB", # MES
                      "PTPRC") # WBC
p2 <- DotPlot(rnaAggr, features = celltype.markers) + theme(axis.text.x = element_text(angle = 90, hjust = 1))


########## Annotation of the clusters for data integration with snATAC dataset ##########
rnaAggr[["orig.clusters"]] <- Idents(object = rnaAggr) # stash cluster idents prior to annotation
new.cluster.ids <- c("PCT","DCT1","CNT","TAL","PCT","TAL","ICA","TAL","PC","PEC","ENDO","PT_KIM1","MES","PODO","DCT2","ICB","ENDO","ENDO")
names(new.cluster.ids) <- levels(rnaAggr)
rnaAggr <- RenameIdents(rnaAggr, new.cluster.ids)
rnaAggr@meta.data$celltype <- rnaAggr@active.ident

p3 <- DimPlot(rnaAggr, reduction = "umap", assay = "SCT", label = TRUE)
p4 <- DotPlot(rnaAggr, features = celltype.markers) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

# create new disease-specific celltype-speciifc identity for FindMarkers function
rnaAggr@meta.data$celltype.stim <- paste(Idents(rnaAggr), rnaAggr@meta.data$disease, sep ="_")
rnaAggr$celltype <- Idents(rnaAggr)
Idents(rnaAggr) <- "celltype.stim"

print("Saving aggregated snRNAseq object as rnaAggr.rds in:")
here("cellrangerRnaAggr")
saveRDS(rnaAggr, file = "rnaAggr.rds")