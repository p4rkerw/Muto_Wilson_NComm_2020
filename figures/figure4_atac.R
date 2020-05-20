library(Seurat)
library(Signac) 
library(monocle3)
library(cicero)
library(BuenColors)
library(EnsDb.Hsapiens.v86)
library(here)
set.seed(1234)

#Use the object with all version of jasper motifs
sub_atac <- readRDS("cellranger_atac_prep/atacAggr_sub97_control_allT.rds")
DefaultAssay(sub_atac) <- "peaks"
count_data <- GetAssayData(sub_atac, slot = "counts")
summ <- summary(count_data)
summ_frame <- data.frame(peak = rownames(count_data)[summ$i],
                         cell.id = colnames(count_data)[summ$j],
                         count = summ$x)

# create cell data set object with cicero constructor
input_cds <- make_atac_cds(summ_frame, binarize = F)
meta <- sub_atac@meta.data
meta$cells <- rownames(meta)
metanames <- rownames(meta)
meta <- merge(input_cds@colData,meta,by.x="cells", by.y="cells")
rownames(meta) <- meta@listData[["cells"]]
input_cds@colData <- meta
input_cds <- detect_genes(input_cds)
input_cds <- estimate_size_factors(input_cds)

input_cds <- preprocess_cds(input_cds, num_dim = 50)
input_cds <- align_cds(input_cds, alignment_group = "orig.ident")
input_cds <- reduce_dimension(input_cds, reduction_method = 'UMAP', 
                              preprocess_method = "Aligned")
plot_cells(input_cds, color_cells_by = "celltype")
input_cds_save <- input_cds #Save the object before removing doublets/low quality cells

input_cds <- cluster_cells(input_cds)
plot_cells(input_cds)
#Remove a couple cluster with various celltypes  (doublets/low quality cells) [12,14,19,21]
input_cds <- input_cds[,clusters(input_cds) %in% c(1:11,13,15:18,20,22,23)]

input_cds <- cluster_cells(input_cds)
input_cds <- learn_graph(input_cds)
input_cds <- order_cells(input_cds) #Choose the three most distant point from PT_KIM1 in PCT and PST for "start points" of the trajectories. 

fig4c_1 <-   plot_cells(input_cds,
                        color_cells_by = "celltype",
                        label_groups_by_cluster=FALSE,
                        label_leaves=FALSE,
                        label_branch_points=FALSE,
                        label_roots=FALSE,
                        label_cell_groups=F,show_trajectory_graph=F,group_label_size=3.6) #png 560x430


cds_subset <- input_cds[,colData(input_cds)$celltype %in% c("PCT","PST","PT_KIM1")] #Subsetting PT and PT_KIM1 clusters 
cds_subset <- choose_cells(cds_subset) #Select the main clucster and exclude the scattered, low-quality cells.

cds_subset <- cluster_cells(cds_subset)
cds_subset <- learn_graph(cds_subset)
cds_subset <- order_cells(cds_subset) #Choose the most distant point from PT_KIM1

fig4c_2 <-   plot_cells(cds_subset,
                      color_cells_by = "pseudotime",
                      label_groups_by_cluster=FALSE,
                      label_leaves=FALSE,
                      label_branch_points=FALSE,
                      label_roots=FALSE) #png 560x430

#marker peaks detection
#markers <- FindMarkers(sub_atac, ident.1 = "PT_KIM1",
#ident.2 = c("PCT","PST"), test.use = 'LR', latent.vars = "nCount_peaks",logfc.threshold = 0.5) # find DN-specific dacs
#cf <- ClosestFeature(rownames(markers), annotation=EnsDb.Hsapiens.v86, sep=c(':','-'))
#markers <- cbind(markers, gene=cf$gene_name, distance=cf$distance)

cds_subset_lin <- cds_subset[,is.finite(pseudotime(cds_subset))]
fig4d_1 <- plot_accessibility_in_pseudotime(cds_subset_lin[c("chr1:100719411-100719996","chr11:26714753-26720418"),],breaks = 8) #VCAM1/SLC5A12
fig4d_2 <- plot_accessibility_in_pseudotime(cds_subset_lin[c("chr15:63040511-63045764","chr4:71338336-71340367"),],breaks = 8) #TPM1/SLC4A4
#png 500x500

p <- pseudotime(cds_subset)
p <- p[order(p)]
p <- as.data.frame(p)
colnames(p) <- "pseudotime"

p$pseudotime_order <- 1:10817
atacPT <- subset(sub_atac, cells = rownames(p)) #Subset PT,PT_KIM1 cells defined in the Cicero object
atacPT <- AddMetaData(atacPT,p)
Idents(atacPT) <- "celltype"

new.cluster.ids <- c("PT","PT","PT_KIM1") #Rename the clusters as just PT vs PT_KIM1 for pseudotime-trajectory.
names(new.cluster.ids) <- levels(atacPT)
atacPT <- RenameIdents(atacPT, new.cluster.ids)
atacPT@meta.data$celltype <- atacPT@active.ident

atacPT <- subset(atacPT,downsample=500) #Down-sizing 500 cells in each celltype
Idents(atacPT) <- "pseudotime_order"
atacPT@active.ident <- ordered(atacPT@active.ident, levels = sort(as.numeric(levels(atacPT))))
motifmatrix <- atacPT@assays[["chromvar"]]@data[rownames(atacPT@assays[["chromvar"]]) %in% c("REL","RELA","FOS::JUN","HNF4A","HNF1A","PPARA::RXRA","NR1H2::RXRA"),]
motifmatrix <- motifmatrix[c("REL","RELA","FOS::JUN","HNF4A","HNF1A","PPARA::RXRA","NR1H2::RXRA"),] #Re-order

#fig4e <- pheatmap::pheatmap(motifmatrix,cluster_cols = F,cluster_rows = F,scale = "column",border_color = "NA",show_colnames=F)

pt <- RunUMAP(object = pt, reduction = 'harmony', dims = 1:40, assay.use = "peaks")
pt <- FindNeighbors(object = pt, reduction = 'harmony', dims = 1:40, assay.use = "peaks")
pt <- FindClusters(object = pt, verbose = FALSE, reduction = 'harmony', assay.use = "peaks")
DefaultAssay(pt) <- "chromvar"
DimPlot(pt, reduction ="umap", label = TRUE, repel = TRUE,group.by = "celltype")+NoLegend()
VlnPlot(pt,"nCount_features",pt.size = 0)
fig4e_1 <- FeaturePlot(sub_atac,features = "RELA",cols =jdb_palette("brewer_yes"))  #png 500x420
fig4e_2 <- FeaturePlot(sub_atac,features = "HNF4A",cols =jdb_palette("brewer_yes")) 
fig4e_3 <- FeaturePlot(pt,features = "RELA",cols =jdb_palette("brewer_yes"))  #png 500x420
fig4e_4 <- FeaturePlot(pt,features = "HNF4A",cols =jdb_palette("brewer_yes")) 

sub_atac <- SetFragments(object = sub_atac, file = "outs_atac/fragments.tsv.gz")

figS11_1 <- CoveragePlot(
  object = sub_atac,
  region = "chr1:100719411-100719996",
  sep = c(":", "-"),
  extend.upstream = 5000,
  extend.downstream = 5000,
  annotation = EnsDb.Hsapiens.v86,
) #VCAM1

figS11_2 <- CoveragePlot(
  object = sub_atac,
  region = "chr15:63040511-63045764",
  sep = c(":", "-"),
  extend.upstream = 5000,
  extend.downstream = 5000,
  annotation = EnsDb.Hsapiens.v86,
) #TPM1

figS11_3 <- CoveragePlot(
  object = sub_atac,
  region = "chr11:26714753-26720418",
  sep = c(":", "-"),
  extend.upstream = 5000,
  extend.downstream = 5000,
  annotation = EnsDb.Hsapiens.v86,
) #SLC5A12

figS11_4 <- CoveragePlot(
  object = sub_atac,
  region = "chr4:71338336-71340367",
  sep = c(":", "-"),
  extend.upstream = 5000,
  extend.downstream = 5000,
  annotation = EnsDb.Hsapiens.v86,
) #SLC4A4

