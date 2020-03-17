library(Seurat)
library(Matrix)
library(openxlsx)
library(dplyr)
library(harmony)
library(DoubletFinder)
library(stringr)
library(tibble)
library(ggplot2)
library(here)


# Doublet removal with the assumption that doublets represent 6% of cells.
## Pre-process Seurat object (standard) --------------------------------------------------------------------------------------
FindDoublets <- function(library_id, seurat_aggregate) {
  obj <- seurat_aggregate
  seurat_obj <- subset(obj, idents = library_id)
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

# read in count matrix downloaded from 
# https://science.sciencemag.org/highwire/filestream/713964/field_highwire_adjunct_files/6/aat1699_DataS1.gz.zip
# archive must be unzipped using a mac utility
counts <- readMM("comparison_human_datasets/young_PMID30093597/tableOfCounts.mtx")
row_labels <- read.table("comparison_human_datasets/young_PMID30093597/tableOfCounts_rowLabels.tsv", header=TRUE)
col_labels <- read.table("comparison_human_datasets/young_PMID30093597/tableOfCounts_colLabels.tsv", header=TRUE)

rownames(counts) <- row_labels$Symbol
colnames(counts) <- col_labels$DropletID

# get barcodes and sampleid that identify control kidney samples
metadata <- read.xlsx("comparison_human_datasets/young_PMID30093597/aat1699-Young-TablesS1-S12-revision1.xlsx",
                      sheet = 'TableS6 - Sample manifest', colNames = TRUE, startRow = 2)
adult_control_kidney_metadata <- dplyr::filter(metadata, TissueDiseaseState == "Normal") %>%
  dplyr::filter(AgeInMonthsPostConception > 60) %>%
  dplyr::filter(Organ == "Kidney") 

adult_control_kidney_sangerid <- dplyr::select(adult_control_kidney_metadata, SangerID)
adult_control_kidney_dropletid <- dplyr::filter(col_labels, SangerID %in% adult_control_kidney_sangerid$SangerID) %>%
  dplyr::select(DropletID)

# filter the count countsrix for barcode id corresponding to control adult kidney samples
counts_control <- counts[,colnames(counts) %in% adult_control_kidney_dropletid$DropletID]

# per the MS the PT-KIM1 / VCAM1 cluster was associated with the cancer samples. Proceed with all samples for seurat analysis.
# create a seurat object with the filtered counts matrix
obj <- CreateSeuratObject(counts_control, min.cells = 3, min.features = 200)
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 20) # manuscript used percent.mt < 20

# take an aggregated snRNA seurat object and a list library_id and find doublets. return a df of doublet barcodes
# send DimPlot and FeaturePlot of doublets for each library to here("plots")
Idents(obj) <- "orig.ident"
list.doublet.bc <- lapply(levels(obj), function(x) {FindDoublets(x, seurat_aggregate = obj)})
doublet_id <- list.doublet.bc %>%
  bind_rows() %>%
  dplyr::rename("doublet_id" = "V2") %>%
  tibble::column_to_rownames(var = "V1") # this is the barcode column
table(doublet_id) # quantify total doublet vs. singlet calls (expect ~5% doublets)
obj <- AddMetaData(obj,doublet_id)

# filter out doublets prior to snRNA preprocessing
Idents(obj) <- "doublet_id"
obj <- subset(obj,idents = "Singlet")


# standard processing pipeline for seurat
obj <- SCTransform(obj, vars.to.regress = c("nCount_RNA"), verbose = TRUE)
obj <- RunPCA(obj)
obj <- RunHarmony(obj, "orig.ident", plot_convergence = TRUE)
obj <- FindNeighbors(obj, reduction = 'harmony')
obj <- FindClusters(obj, reduction = 'harmony')
obj <- RunUMAP(obj, reduction = 'harmony', dims = 1:30)

p1 <- DimPlot(obj, label = T, reduction = "umap") + NoLegend()



celltype.markers <- c("CUBN","LRP2","HAVCR1","VCAM1", # PT and PT-KIM1+ markers
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
p2 <- DotPlot(obj, features = celltype.markers) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
CombinePlots(list(p1,p2))

# perform label transfer from reference human adult kidney dataset
rnaAggr <- readRDS("cellranger_rna_prep/rnaAggr_control.rds")
kidney.anchors <- FindTransferAnchors(reference = rnaAggr, query = obj, 
                                        dims = 1:30)
predictions <- TransferData(anchorset = kidney.anchors, refdata = rnaAggr$celltype, 
                            dims = 1:30)
obj <- AddMetaData(obj, metadata = predictions)

Idents(obj) <- "predicted.id"
DimPlot(obj, label = TRUE, reduction="umap")
