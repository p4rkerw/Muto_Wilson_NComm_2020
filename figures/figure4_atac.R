library(Seurat)
library(Signac) 
library(monocle3)
library(cicero)
library(BuenColors)
library(EnsDb.Hsapiens.v86)
library(here)
set.seed(1234)

#Use the object with all version of jasper motifs
sub_atac <- readRDS("cellranger_atac_prep/sub_atac_sub97_control_allT.rds")
atacPT <- subset(sub_atac,idents = c("PCT","PST","PT_KIM1"))
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
plot_cells(input_cds, color_cells_by = "orig.ident")
plot_cells(input_cds, color_cells_by = "celltype")
plot_cells(input_cds, color_cells_by = "predicted.id")
input_cds_save <- input_cds
#Remove a couple cluster with various celltypes
input_cds <- cluster_cells(input_cds)
plot_cells(input_cds2)
#Remove a couple cluster with various celltypes (12,14,19,21)
input_cds2 <- input_cds[,clusters(input_cds) %in% c(1:11,13,15:18,20,22,23)]

input_cds <- cluster_cells(input_cds)
input_cds <- learn_graph(input_cds,close_loop = T)
input_cds <- order_cells(input_cds) #Choose the most distant point from PT_KIM1
plot_cells(input_cds, color_cells_by = "orig.ident",cell_size = 1,label_cell_groups=F)
plot_cells(input_cds, color_cells_by = "celltype",cell_size = 1,label_cell_groups=F)
plot_cells(input_cds, color_cells_by = "pseudotime",cell_size = 1,trajectory_graph_color="yellow")
plot_cells(input_cds, color_cells_by = "prediction.score.PT",cell_size = 1,label_cell_groups=T)

cds_subset <- input_cds[,colData(input_cds)$celltype %in% c("PCT","PST","PT_KIM1")]
cds_subset <- choose_cells(cds_subset) #Select the main clucster

cds_subset <- cluster_cells(cds_subset)
cds_subset <- learn_graph(cds_subset,close_loop = T)
cds_subset <- order_cells(cds_subset) #Choose the most distant point from PT_KIM1

fig4e_1 <-   plot_cells(cds_subset,
                        color_cells_by = "celltype",
                        label_groups_by_cluster=FALSE,
                        label_leaves=FALSE,
                        label_branch_points=FALSE,
                        label_roots=FALSE,
                        ,label_cell_groups=F,show_trajectory_graph=F,group_label_size=3.6) #png 560x430

fig4e_2 <-   plot_cells(cds_subset,
                      color_cells_by = "pseudotime",
                      label_groups_by_cluster=FALSE,
                      label_leaves=FALSE,
                      label_branch_points=FALSE,
                      label_roots=FALSE) #png 560x430

#marker peaks detection
markers <- FindMarkers(sub_atac, ident.1 = "PT_KIM1",
ident.2 = c("PCT","PST"), test.use = 'LR', latent.vars = "nCount_peaks",logfc.threshold = 0.5) # find DN-specific dacs
cf <- ClosestFeature(rownames(markers), annotation=EnsDb.Hsapiens.v86, sep=c(':','-'))
markers <- cbind(markers, gene=cf$gene_name, distance=cf$distance)

cds_subset_lin <- cds_subset[,is.finite(pseudotime(cds_subset))]
fig4f_1 <- plot_accessibility_in_pseudotime(cds_subset_lin[c("chr1:100719411-100719996","chr15:63040511-63045764"),],breaks = 8) #VCAM1/TNIK
fig4f_2 <- plot_accessibility_in_pseudotime(cds_subset_lin[c("chr11:26714753-26720418","chr4:71338336-71340367"),],breaks = 8) #SLC5A12/SLC4A4

#png 467x639
p <- pseudotime(cds_subset)
p <- p[order(p)]
p <- as.data.frame(p)
colnames(p) <- "pseudotime"
#p$pseudotime_order <- c(rep(1:100, each=10)) #8264 cells
#p$pseudotime_order  <- cut(p$pseudotime, breaks = seq(-0.0001,8.303220,length = 11),label=1:10)
p$pseudotime_order <- 1:10817
atacPT <- subset(sub_atac, cells = rownames(p))
atacPT <- AddMetaData(atacPT,p)
Idents(atacPT) <- "celltype"

new.cluster.ids <- c("PT","PT","PT_KIM1")
names(new.cluster.ids) <- levels(atacPT)
atacPT <- RenameIdents(atacPT, new.cluster.ids)
atacPT@meta.data$celltype <- atacPT@active.ident
atacPT2 <- subset(atacPT,downsample=500)
Idents(atacPT2) <- "pseudotime_order"
atacPT2@active.ident <- ordered(atacPT2@active.ident, levels = sort(as.numeric(levels(atacPT2))))
#peaks_average <- AverageExpression(atacPT,assays = "chromvar")
#aver_chromvar  <- AverageExpression(atacPT2,assays = "chromvar",features = c("REL","RELA","FOS::JUN","HNF4A","HNF1A","PPARA::RXRA","NR1H2::RXRA"))
#pheatmap::pheatmap(peaks_average[["chromvar"]],cluster_cols = F,cluster_rows = F,scale = "column",border_color = "NA") 
#DimPlot(atacPT)
motifmatrix <- atacPT2@assays[["chromvar"]]@data[rownames(atacPT2@assays[["chromvar"]]) %in% c("REL","RELA","FOS::JUN","HNF4A","HNF1A","PPARA::RXRA","NR1H2::RXRA"),]
motifmatrix <- motifmatrix[c(1,2,4,5,6,3,7),]
fig4g <- pheatmap::pheatmap(motifmatrix,cluster_cols = F,cluster_rows = F,scale = "column",border_color = "NA",show_colnames=F)

FeaturePlot(atacPT,features = "HNF1A")
