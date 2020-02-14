# this script will eliminate doublets from an aggregated snRNA object prior to preprocessing
library(Seurat) # 3.02
library(ggplot2)
library(harmony) # 1.0
library(Rcpp)
library(DoubletFinder)
library(openxlsx)
library(here)
library(dplyr)
library(stringr)
library(tibble)
set.seed(1234)

outs <- "cellranger_rna_aggr_control/outs/"

# load aggregated snRNAseq data from cellranger aggregate matrix and create seurat objects
counts <- Read10X(here(outs,"filtered_feature_bc_matrix"))
metadata <- read.csv(here(outs,"aggregation.csv"))

# load aggregated snRNAseq data from cellranger aggregate matrix and create seurat objects
rnaAggr <- CreateSeuratObject(counts = counts, min.cells = 10, min.features = 500, 
                              project = "RNA")
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

VlnPlot(object = rnaAggr, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "orig.ident",pt.size=0)
VlnPlot(object = rnaAggr, features = c("percent.rps", "percent.rpl"), ncol = 2, group.by = "orig.ident",pt.size=0)

# filter the aggregated dataset for low quality cells
rnaAggr <- subset(rnaAggr, subset = nFeature_RNA > 500 # use same filter params as aggr
                  & nFeature_RNA < 4000 
                  & nCount_RNA < 16000 
                  & percent.mt < 0.8 
                  & percent.rps < 0.4 
                  & percent.rpl < 0.4)

# Doublet removal with the assumption that doublets represent 6% of cells.
## Pre-process Seurat object (standard) --------------------------------------------------------------------------------------

FindDoublets <- function(library_id, seurat_aggregate) {
  rnaAggr <- seurat_aggregate
  seurat_obj <- subset(rnaAggr, idents = library_id)
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- ScaleData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
  seurat_obj <- RunPCA(seurat_obj)
  # ElbowPlot(seurat_obj)
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)
  DimPlot(seurat_obj)

  ## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
  sweep.res.list_kidney <- paramSweep_v3(seurat_obj, PCs = 1:20, sct = F)
  sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = FALSE)
  bcmvn_kidney <- find.pK(sweep.stats_kidney)
  pK <- bcmvn_kidney %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
    filter(BCmetric == max(BCmetric)) %>%
    select(pK) 
  pK <- as.numeric(as.character(pK[[1]]))
  seurat_doublets <- doubletFinder_v3(seurat_obj, PCs = 1:20, pN = 0.25, pK = pK,
                               nExp = round(0.05*length(seurat_obj@active.ident)), 
                               reuse.pANN = FALSE, sct = F)
 
  # create doublet groupings and visualize results
  DF.class <- names(seurat_doublets@meta.data) %>% str_subset("DF.classifications")
  pANN <- names(seurat_doublets@meta.data) %>% str_subset("pANN")
  
  p1 <- ggplot(bcmvn_kidney, aes(x=pK, y=BCmetric)) +
    geom_bar(stat = "identity") + 
    ggtitle(paste0("pKmax=",pK)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  p2 <- DimPlot(seurat_doublets, group.by = DF.class)
  p3 <- FeaturePlot(seurat_doublets, features = pANN)
  
  outFile <- paste0("doublets.",library_id,".pdf")
  dir.create("plots", showWarnings = FALSE)
  pdf(here("plots",outFile))
    print(p1) # need to use print() when drawing pdf in a function call
    print(p2)
    print(p3)
  dev.off()
  
  # create a df of barcodes and doublet designations
  df_doublet_barcodes <- as.data.frame(cbind(rownames(seurat_doublets@meta.data), seurat_doublets@meta.data[[DF.class]]))
  return(df_doublet_barcodes)
}


# take an aggregated snRNA seurat object and a list library_id and find doublets. return a df of doublet barcodes
# send DimPlot and FeaturePlot of doublets for each library to here("plots")
Idents(rnaAggr) <- "orig.ident"
list.doublet.bc <- lapply(orig.ident, function(x) {FindDoublets(x, seurat_aggregate = rnaAggr)})
doublet_id <- list.doublet.bc %>%
  bind_rows() %>%
  dplyr::rename("doublet_id" = "V2") %>%
  tibble::column_to_rownames(var = "V1") # this is the barcode column
table(doublet_id) # quantify total doublet vs. singlet calls (expect ~5% doublets)
  
# add doublet calls to aggregated snRNA object as doublet_id in meta.data slot
rnaAggr <- AddMetaData(rnaAggr,doublet_id)
# remove(doublet.meta_data)

# filter out doublets prior to snRNA preprocessing
Idents(rnaAggr) <- "doublet_id"
rnaAggr <- subset(rnaAggr,idents = "Singlet")
saveRDS(rnaAggr,here("cellranger_rna_prep","rnaAggr_control_nodoublets.rds"))


# run sctransform
# Regress out the mitochondrial reads and nCount_RNA
rnaAggr <- SCTransform(rnaAggr, vars.to.regress = c("percent.mt","percent.rpl", "percent.rps","nCount_RNA"), verbose = TRUE)
rnaAggr <- RunPCA(rnaAggr, verbose = TRUE)
# ElbowPlot(rnaAggr, ndims = 50) # to determine number of dimensions for clustering
rnaAggr <- RunHarmony(rnaAggr, "orig.ident", plot_convergence = TRUE)
rnaAggr <- FindNeighbors(rnaAggr, dims = 1:24, verbose = TRUE, reduction = "harmony")
rnaAggr <- FindClusters(rnaAggr, verbose = TRUE, resolution = 0.7, reduction = "harmony")
rnaAggr <- RunUMAP(rnaAggr, dims = 1:24, verbose = TRUE, reduction = "harmony")

# visualize the clustering
# DimPlot(rnaAggr, reduction = "UMAP", assay = "SCT")
p1 <- DimPlot(rnaAggr, reduction = "umap", assay = "SCT", label = TRUE) + ggtitle("snRNA Seurat Clustering with Harmony No Doublets")

celltype.markers <- c("CUBN","HAVCR1","SLC5A1","SLC5A2", # PT and PT-KIM1+ markers
                      "CFH", # PEC
                      "SLC12A1", # TAL NKCC2
                      "SLC12A3","TRPM6", # DCT1 and DCT2 NCC
                      "SCNN1G","TRPV5", # DCT2/CNT ENaC
                      "CALB1", # CNT
                      "AQP2", # PC
                      "ATP6V0D2", # ICA and ICB
                      "SLC4A1","SLC26A7", # ICA
                      "SLC26A4", # ICB
                      "NPHS1","NPHS2", # PODO
                      "PECAM1","FLT1","EMCN", # ENDO
                      "CLDN5", # GEC
                      "ITGA8","PDGFRB", # MES
                      "ACTA2","CALD1", # FIB
                      "PTPRC") # WBC
p2 <- DotPlot(rnaAggr, features = celltype.markers) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

CombinePlots(plots = list(p1, p2))
########## Annotation of the clusters for data integration with snATAC dataset ##########
rnaAggr[["orig.clusters"]] <- Idents(object = rnaAggr) # stash cluster idents prior to annotation
new.cluster.ids <- c("PCT","DCT1","TAL","CNT","PCT",
                     "TAL","ICA","TAL","PC","ENDO",
                     "PEC","DCT2","PODO","PT_KIM1","ICB",
                     "ENDO","MES","FIB","ENDO","LEUK")
names(new.cluster.ids) <- levels(rnaAggr)
rnaAggr <- RenameIdents(rnaAggr, new.cluster.ids)

# reorder the idents and save celltype annotations in celltype slot metadata
levels(rnaAggr) <- c("PCT","PT_KIM1","PEC","TAL","DCT1",
                     "DCT2","CNT","PC","ICA","ICB",
                     "PODO","ENDO","MES","FIB","LEUK")
rnaAggr@meta.data$celltype <- rnaAggr@active.ident

# create low-resolution celltype identities for snATAC thresholding (ie group PT and PT-KIM1 and distal nephron together)
lowres.cluster.ids <- c("PCT","PCT","PEC","TAL","DCT_CNT_PC",
                        "DCT_CNT_PC","DCT_CNT_PC","DCT_CNT_PC","ICA","ICB",
                        "PODO","ENDO","MES_FIB","MES_FIB","LEUK")
names(lowres.cluster.ids) <- levels(rnaAggr)
rnaAggr <- RenameIdents(rnaAggr, lowres.cluster.ids)
rnaAggr[["lowres.celltype"]] <- Idents(rnaAggr)

# reset celltype as primary ident before saving and plotting
Idents(rnaAggr) <- "celltype"

# redraw umap and dotplot with reordered idents
p3 <- DimPlot(rnaAggr, reduction = "umap", assay = "SCT", label = TRUE) + ggtitle("snRNA Annotated Celltypes")
p4 <- DotPlot(rnaAggr, features = celltype.markers) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

CombinePlots(list(p3,p4))
# draw pdf plots for before and after annotation
dir.create("plots", showWarnings = FALSE)
pdf(here("plots","umap.rnaAggr.pdf"))
p1
p2
p3
p4
dev.off()

# save preprocessed rnAggr file
print("Saving aggregated snRNAseq object as rnaAggr.rds in:")
dir.create("cellranger_rna_prep", showWarnings = FALSE)
here("cellranger_rna_prep")
saveRDS(rnaAggr, file = here("cellranger_rna_prep","rnaAggr_control.rds"))

        
        
        
        
        
        
        
        
        
        
        
        