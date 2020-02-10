# this script will preprocess aggregated snATACseq data from 5 healthy control and 5 diabetic kidney cortex samples
# counted and aggregated by cellranger-atac v1.2.0 without library normalization

library(Signac) #version 0.1.5
library(Seurat) #version 3.1.1
library(GenomeInfoDb)
library(harmony) #version 1.0
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(here)
library(tibble)
library(dplyr)
set.seed(1234)


# initialize output directories for aggregated and individual snATACseq data
atacAggrDir <- "cr12AtacAggr10/outs"
atacCountDir <- here()
rnaAggrDir <- "cellrangerRnaAggr"
fragment.path <- here("cr12AtacAggr10/outs", "fragments.tsv.gz")

# load aggregated snATACseq data obtained from cellrangerAtacAggr.sh v1.2.0 and create a seurat object
counts <- Read10X_h5(here(atacAggrDir,"filtered_peak_bc_matrix.h5"))
metadata <- read.csv(here(atacAggrDir, "singlecell.csv"), header = TRUE, row.names = 1)
aggcsv <- read.csv(here(atacAggrDir, "aggregation_csv.csv"),header = TRUE, row.names = 1)
atacAggr <- CreateSeuratObject(counts = counts, assay = 'peaks', project = 'ATAC',
                               min.cells = 5, meta.data = metadata)
atacAggr <- SetFragments(object = atacAggr, file = fragment.path)
atacAggr
rm(counts, metadata)

# Add the patient information and disease status to the metadata of the Seurat object
gemgroup <- sapply(strsplit(rownames(atacAggr@meta.data), split="-"), "[[", 2) 
current.gemgroups <- seq(length(rownames(aggcsv))) # no. gemgroups is no. samples
orig.ident <- rownames(aggcsv)
sampleID <- plyr::mapvalues(gemgroup, from = current.gemgroups, to = orig.ident)
atacAggr <- AddMetaData(object=atacAggr, metadata=data.frame(orig.ident=sampleID, 
  row.names=rownames(atacAggr@meta.data)))

# Create an additional metadata grouping for control vs diabetes in meta.data$disease
disease.groups <- aggcsv$group # control vs. diabetes
groupID <-plyr::mapvalues(atacAggr@meta.data$orig.ident,
  from = orig.ident,to=disease.groups)
levels(groupID) <- levels(disease.groups)
atacAggr <- AddMetaData(atacAggr, metadata=data.frame(disease=groupID,
    row.names=rownames(atacAggr@meta.data)))


# calculating nucleosome signal
atacAggr <- NucleosomeSignal(object = atacAggr) #Preparation of nucleosome_signal

# Add the metadata for mitochondrial fragments from individual snATACseq counts into the aggregated dataset 
# Collect total fragment number data from each of the original CellRangerATAC datasets.
metaqc <- 
  lapply(current.gemgroups, 
    function(gemgroup) {
      sampleName <- rownames(aggcsv)[gemgroup]
      df <- read.csv(file = paste0(sampleName, "/outs/singlecell.csv"), header=TRUE, row.names=1)
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
atacAggr <- RunTFIDF(atacAggr)
atacAggr <- FindTopFeatures(atacAggr, min.cutoff = 'q1')
atacAggr <- RunSVD(
  object = atacAggr,
  assay = 'peaks',
  reduction.key = 'PCA_',
  reduction.name = 'pca'
)
atacAggr <- RunUMAP(object = atacAggr, reduction = "pca", dims = 1:12)
atacAggr <- FindNeighbors(object = atacAggr, reduction = "pca", dims = 1:12)
atacAggr <- FindClusters(object = atacAggr, verbose = FALSE)
atacAggr <- RunHarmony(atacAggr, group.by.vars = "orig.ident", assay.use = "peaks")
atacAggr <- RunUMAP(object = atacAggr, reduction = "harmony", dims = 1:33,
                    assay.use = "peaks")

# saveRDS(atacAggr, file = "atacAggr.rds")

# create a gene activity matrix to estimate transcriptional activity from snATACseq data
# extract gene coordinates from Ensembl, and ensure name formatting is consistent with Seurat object 
gene.coords <- genes(EnsDb.Hsapiens.v86, filter = ~ gene_biotype == "protein_coding")
seqlevelsStyle(gene.coords) <- 'UCSC'
genebody.coords <- keepStandardChromosomes(gene.coords, pruning.mode = 'coarse')
genebodyandpromoter.coords <- Extend(x = gene.coords, upstream = 2000, downstream = 0)

# Memory-exhaustive step:
gene.activities <- FeatureMatrix(fragments = fragment.path, features = genebodyandpromoter.coords, cells = colnames(atacAggr), chunk = 10)
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

# saveRDS(atacAggr, file = "atacAggr.rds")

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
rnaAggr <- readRDS(here(rnaAggrDir,"rnaAggr.rds"))

# identify anchors to transfer cell labels from snRNAseq to snATACseq "RNA" gene activity scores
transfer.anchors <- FindTransferAnchors(
  reference = rnaAggr,
  query = atacAggr,
  reference.assay = 'RNA',
  query.assay = 'RNA',
  reduction = 'cca',
)

rnaAggr@active.ident -> rnaAggr@meta.data$celltype
predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = rnaAggr$celltype,
  weight.reduction = atacAggr[["lsi"]]
)

# add predicted cell types to the snATACseq object
atacAggr <- AddMetaData(atacAggr, metadata = predicted.labels)
atacAggr$predicted.id <- factor(atacAggr$predicted.id, levels = levels(rnaAggr))

p1 <- DimPlot(atacAggr, group.by = "predicted.id", label = TRUE, repel = TRUE) + ggtitle("snATAC-seq cells") + 
    NoLegend() + scale_colour_hue(drop = FALSE)
p2 <- DimPlot(rnaAggr, group.by = "celltype", label = TRUE, repel = TRUE) + ggtitle("snRNA-seq cells") + 
    NoLegend()
CombinePlots(plots = list(p1, p2))

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

p1 <- DimPlot(atacAggr, group.by = "predicted.id", label = TRUE, repel = TRUE) + ggtitle("snATAC-seq cells") + 
    NoLegend() + scale_colour_hue(drop = FALSE)
p2 <- DimPlot(rnaAggr, label = TRUE, repel = TRUE) + ggtitle("snRNA-seq cells") + 
    NoLegend()
CombinePlots(plots = list(p1, p2))

print("Saving integrated snRNAseq - snATACseq file to:")
getwd()
saveRDS(atacAggr, file = "atacAggr.rds")