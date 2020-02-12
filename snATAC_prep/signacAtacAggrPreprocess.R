# this script will preprocess aggregated snATACseq data from 5 healthy control kidney cortex samples
# counted and aggregated by cellranger-atac v1.2.0 without library normalization

library(Signac) #version 0.1.5
library(Seurat) #version 3.0.2
library(GenomeInfoDb)
library(harmony) #version 1.0
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(here)
library(tibble)
library(dplyr)
set.seed(1234)


# initialize output directories for aggregated and individual snATACseq data
outs_atac <- "cellranger_atac_aggr_control/outs"
prep_rna <- "cellranger_rna_prep"

# load aggregated snATACseq data obtained from cellrangerAtacAggr.sh v1.2.0 and create a seurat object
counts <- Read10X_h5(here(outs_atac,"filtered_peak_bc_matrix.h5"))
metadata <- read.csv(here(outs_atac, "singlecell.csv"), header = TRUE, row.names = 1)
aggcsv <- read.csv(here(outs_atac, "aggregation_csv.csv"),header = TRUE, row.names = 1)
atacAggr <- CreateSeuratObject(counts = counts, assay = 'peaks', project = 'ATAC',
                               min.cells = 5, meta.data = metadata)
atacAggr <- SetFragments(object = atacAggr, file = here(outs_atac, "fragments.tsv.gz"))
atacAggr
rm(counts, metadata)

# Add the patient information and disease status to the metadata of the Seurat object
gemgroup <- sapply(strsplit(rownames(atacAggr@meta.data), split="-"), "[[", 2) 
current.gemgroups <- seq(length(rownames(aggcsv))) # no. gemgroups is no. samples
orig.ident <- rownames(aggcsv)
sampleID <- plyr::mapvalues(gemgroup, from = current.gemgroups, to = orig.ident)
atacAggr <- AddMetaData(object=atacAggr, metadata=data.frame(orig.ident=sampleID, 
  row.names=rownames(atacAggr@meta.data)))

# # Create an additional metadata grouping for control vs diabetes in meta.data$disease
# disease.groups <- aggcsv$group # control vs. diabetes
# groupID <-plyr::mapvalues(atacAggr@meta.data$orig.ident,
#   from = orig.ident,to=disease.groups)
# levels(groupID) <- levels(disease.groups)
# atacAggr <- AddMetaData(atacAggr, metadata=data.frame(disease=groupID,
#     row.names=rownames(atacAggr@meta.data)))


# calculating nucleosome signal
atacAggr <- NucleosomeSignal(object = atacAggr)

# Add the metadata for mitochondrial fragments from individual snATACseq counts into the aggregated dataset 
# Collect total fragment number data from each of the original CellRangerATAC datasets.
metaqc <- 
  lapply(current.gemgroups, 
    function(gemgroup) {
      sampleName <- rownames(aggcsv)[gemgroup]
      file <- here("cellranger_atac_counts","version_1.2",sampleName,"outs","singlecell.csv")
      df <- read.csv(file, header=TRUE, row.names=1)
      rownames(df) <- paste(substr(rownames(df), 1, 16), gemgroup, sep = "-") # change gemgroup to reflect sample order
      df <- tibble::rownames_to_column(df, var = "barcode")
      return(df)
    }
  ) %>%
  bind_rows() %>%
  tibble::column_to_rownames(var = "barcode") %>%
  dplyr::select(c("total","mitochondrial"))
atacAggr <- AddMetaData(atacAggr,metaqc)
remove(metaqc)

# QC metrics and filtering
# atacAggr <- subset(atacAggr , subset = nucleosome_signal < 50) # remove outliers
atacAggr$pct_reads_in_peaks <- atacAggr$peak_region_fragments / atacAggr$passed_filters * 100 # %peaks/fragments
atacAggr$blacklist_ratio <- atacAggr$blacklist_region_fragments / atacAggr$peak_region_fragments # %blacklistregion/peaks
atacAggr$mito_ratio <- atacAggr$mitochondrial / atacAggr$total # %peaks/fragments * 100
atacAggr$nucleosome_group <- ifelse(atacAggr$nucleosome_signal > 10, 'NS > 10', 'NS < 10')

# filter the aggregated snATACseq object using empirically-determined QC parameters
atacAggr <- subset(atacAggr, subset = peak_region_fragments > 2500 & peak_region_fragments < 25000 &
 pct_reads_in_peaks > 15 & blacklist_ratio < 0.001 & nucleosome_signal < 4 & mito_ratio < 0.25)

# perform normalization and dimensional reduction of aggregated snATACseq object
# call the singular value decomposition reduction "PCA" so harmony can recognize it
# this is actually an LSI reduction
atacAggr <- RunTFIDF(atacAggr)
atacAggr <- FindTopFeatures(atacAggr, min.cutoff = 'q1')
atacAggr <- RunSVD(
  object = atacAggr,
  assay = 'peaks',
  reduction.key = 'PCA_', # reduction is named "PCA" but is actually LSI
  reduction.name = 'pca'
)
# ElbowPlot(atacAggr, ndim = 40) # select number of dimensions for UMAP embedding
atacAggr <- RunUMAP(object = atacAggr, reduction = "pca", dims = 1:12)
atacAggr <- FindNeighbors(object = atacAggr, reduction = "pca", dims = 1:12)
atacAggr <- FindClusters(object = atacAggr, verbose = TRUE)
atacAggr <- RunHarmony(atacAggr, group.by.vars = "orig.ident", assay.use = "peaks")
atacAggr <- RunUMAP(object = atacAggr, reduction = "harmony", dims = 1:20,
                    assay.use = "peaks")

# plot the original clustering prior to label transfer
p1 <- DimPlot(atacAggr) + ggtitle("Original ATAC Clustering")

dir.create("cellranger_atac_prep")
# saveRDS(atacAggr, file = here("cellranger_atac_prep", "atacAggr_control.rds"))

# create a gene activity matrix to estimate transcriptional activity from snATACseq data
# extract gene coordinates from Ensembl, and ensure name formatting is consistent with Seurat object 
gene.coords <- genes(EnsDb.Hsapiens.v86, filter = ~ gene_biotype == "protein_coding")
seqlevelsStyle(gene.coords) <- 'UCSC'
genebody.coords <- keepStandardChromosomes(gene.coords, pruning.mode = 'coarse')
genebodyandpromoter.coords <- Extend(x = gene.coords, upstream = 2000, downstream = 0)

# Memory-intensive step:
gene.activities <- FeatureMatrix(fragments = here(outs_atac, "fragments.tsv.gz"), features = genebodyandpromoter.coords, cells = colnames(atacAggr), chunk = 10)
gene.key <- genebodyandpromoter.coords$gene_name
names(gene.key) <- GRangesToString(grange = genebodyandpromoter.coords)
rownames(gene.activities) <- gene.key[rownames(gene.activities)]

# add the gene activity matrix to the Seurat snATACseq object as a new assay, and normalize it
atacAggr[['RNA']] <- CreateAssayObject(counts = gene.activities)
atacAggr <- NormalizeData(
  object = atacAggr,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(atacAggr$nCount_RNA)
)

# saveRDS(atacAggr, file = here("cellranger_atac_prep", "atacAggr_control.rds"))

# -------------- START LABEL TRANSFER ---------------------
# create an RNA assay based on gene activity scores, process, and embed with UMAP
DefaultAssay(atacAggr) <- 'RNA'
atacAggr <- RunTFIDF(atacAggr)
atacAggr <- FindTopFeatures(atacAggr, min.cutoff = 'q1')
atacAggr <- RunSVD(
  object = atacAggr,
  assay = 'peaks',
  reduction.key = 'LSI_',
  reduction.name = 'lsi'
)
atacAggr  <- RunUMAP(object = atacAggr, reduction = "lsi", dims = 1:12)

# load the corresponding aggregated snRNAseq object
rnaAggr <- readRDS(here(prep_rna,"rnaAggr_control.rds"))

# identify anchors to transfer cell labels from snRNAseq to snATACseq "RNA" gene activity scores
transfer.anchors <- FindTransferAnchors(
  reference = rnaAggr,
  query = atacAggr,
  reference.assay = 'RNA',
  query.assay = 'RNA',
  reduction = 'cca',
)

# this is high-resolution celltype prediction which is great for predicting celltypes
# but may is not best for thresholding snATAC data
rnaAggr@active.ident -> rnaAggr@meta.data$celltype
predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = rnaAggr$celltype,
  weight.reduction = atacAggr[["lsi"]]
)

# add predicted cell types to the snATACseq object
atacAggr <- AddMetaData(atacAggr, metadata = predicted.labels)
atacAggr$predicted.id <- factor(atacAggr$predicted.id, levels = levels(rnaAggr@meta.data$celltype))
atacAggr$highres.predicted.id <- atacAggr$predicted.id # save the predicted celltype o/w predicted.id is overwritten

# this is low-resolution celltype prediction which may be better for thresholding snATAC data
lowres.predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = rnaAggr$lowres.celltype,
  weight.reduction = atacAggr[["lsi"]]
)

atacAggr <- AddMetaData(atacAggr, metadata = lowres.predicted.labels)
atacAggr$lowres.predicted.id <- factor(atacAggr$predicted.id, levels = levels(rnaAggr@meta.data$lowres.celltype))


p2 <- DimPlot(atacAggr, group.by = "lowres.predicted.id", label = TRUE, repel = TRUE) + ggtitle("snATAC-seq cells Before Harmony") + 
    NoLegend() + scale_colour_hue(drop = FALSE)
p3 <- DimPlot(rnaAggr, reduction = "umap", assay = "SCT", label = TRUE, repel = TRUE) + ggtitle("snRNA-seq Annotated Celltypes") + 
    NoLegend()
CombinePlots(plots = list(p2, p3))

# perform batch effect correction with harmony
atacAggr <- RunSVD(
  object = atacAggr,
  assay = 'peaks',
  reduction.key = 'pca_',
  reduction.name = 'pca')
atacAggr <- RunHarmony(atacAggr, "orig.ident", plot_convergence = TRUE, assay.use = "peaks")
atacAggr <- RunUMAP(object = atacAggr, reduction = 'harmony', dims = 1:33, assay.use = "peaks")
atacAggr <- FindNeighbors(object = atacAggr, reduction = 'harmony', dims = 1:33, assay.use = "peaks")
atacAggr <- FindClusters(object = atacAggr, verbose = FALSE, reduction = 'harmony', assay.use = "peaks")

p4 <- DimPlot(atacAggr, group.by = "highres.predicted.id", label = TRUE, repel = TRUE) + ggtitle("snATAC-seq Predicted Celltypes After Harmony") + 
    NoLegend() + scale_colour_hue(drop = FALSE)
p5 <- DimPlot(rnaAggr, reduction = "umap", assay = "SCT", label = TRUE, repel = TRUE) + ggtitle("snRNA-seq Annotated Celltypes") + 
  NoLegend()
CombinePlots(plots = list(p4, p5))

# plot snATAC data with increasingly stringent predicted.id thresholds
p6 <- hist(atacAggr@meta.data$prediction.score.max, main = "Prediction Score for snATAC")
sub_atac <- subset(atacAggr, subset = prediction.score.max > 0.95)

p7 <- DimPlot(sub_atac, group.by = "lowres.predicted.id", label = TRUE, repel = TRUE) + ggtitle("snATAC-seq Predicted Celltypes After Harmony \n and 95% Prediction Threshold Using \n Low-res Celltypes") + 
  NoLegend() + scale_colour_hue(drop = FALSE)

CombinePlots(plots = list(p4, p7))

# remove doublet barcodes identified by doubletfinder and redo batch correction and clustering

print("Saving integrated snRNAseq - snATACseq file to:")
suppressWarnings(dir.create("cellranger_atac_prep"))
here("cellranger_atac_prep")
saveRDS(atacAggr, file = here("cellranger_atac_prep", "atacAggr_control.rds"))
