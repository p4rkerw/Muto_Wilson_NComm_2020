library(Seurat) # 3.0.2
library(monocle3)
library(ggplot2)
library(dplyr)
library(Matrix)
library(here)
set.seed(1234)

#snRNA-seq data
rnaAggr <- readRDS("cellranger_rna_prep/rnaAggr_control.rds")
#Idents(rnaAggr) <- "celltype"
count_data <- GetAssayData(rnaAggr, assay = "RNA",slot = "counts")

gene_metadata <- as.data.frame(rownames(count_data))
colnames(gene_metadata) <- "gene_short_name"
rownames(gene_metadata) <- gene_metadata$gene_short_name

cds <- new_cell_data_set(as(count_data, "sparseMatrix"),
                         cell_metadata = rnaAggr@meta.data,
                         gene_metadata = gene_metadata)

cds <- preprocess_cds(cds, num_dim = 100)
cds = align_cds(cds, num_dim = 100, alignment_group = "orig.ident")
cds = reduce_dimension(cds,preprocess_method = "Aligned")
fig4a <- plot_cells(cds, color_cells_by="celltype", 
           group_label_size = 3)
cds_subset <- choose_cells(cds) #subseting PT for later analysis: before adding pseudotime

#whole dataset
cds <- cluster_cells(cds)
cds <- learn_graph(cds)
plot_cells(cds,
           color_cells_by = "celltype",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

cds <- order_cells(cds)

fig4b <-   plot_cells(cds,
           color_cells_by = "pseudotime",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           label_roots=FALSE) #png 560x430

plot_cells(cds_subset,
           color_cells_by = "celltype",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,label_cell_groups=T,show_trajectory_graph=F,group_label_size=3.6) #png 460x430

#subclustering

cds_subset <- choose_cells(cds)
cds_subset <- cds_subset[,colData(cds_subset)$celltype %in% c("PT","PT_KIM1")]
cds_subset <- cluster_cells(cds_subset)
cds_subset <- learn_graph(cds_subset,close_loop = FALSE)
cds_subset <- order_cells(cds_subset)
fig4c_1 <- plot_cells(cds_subset,
           color_cells_by = "pseudotime",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           show_trajectory_graph=T) #png 560x430
fig4c_2 <-  plot_cells(cds_subset,
           genes="HAVCR1",cell_size = 1,
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)

genes <- c("VCAM1","ITGB8","TPM1","LRP2","CUBN","SLC5A12")
lineage_cds <- cds_subset[rowData(cds_subset)$gene_short_name %in% genes,]
fig4d <- plot_genes_in_pseudotime(lineage_cds,
                         color_cells_by="celltype",
                         min_expr=1,
                         cell_size=0.1,
                         trend_formula = "~ splines::ns(pseudotime, df=10)",
                         panel_order = c("VCAM1","ITGB8","TPM1","LRP2","CUBN","SLC5A12")
                         ) #png 500x670



# create cell data set object with cicero constructor
input_cds <- make_atac_cds(summ_frame, binarize = T)
meta <- atacPT@meta.data
meta$cells <- rownames(meta)
metanames <- rownames(meta)
meta <- merge(input_cds@colData,meta,by.x="cells", by.y="cells")
rownames(meta) <- meta@listData[["cells"]]
input_cds@colData <- meta
input_cds <- detect_genes(input_cds)
input_cds <- estimate_size_factors(input_cds)
#marker peaks detection
Idents(atacPT) <- "celltype"
markers <- FindMarkers(atacPT, ident.1 = "PT(KIM1+)",
                       ident.2 = "PT", test.use = 'LR', latent.vars = "nCount_peaks") # find DN-specific dacs
cf <- ClosestFeature(rownames(markers), annotation=EnsDb.Hsapiens.v86, sep=c(':','-'))
markers <- cbind(markers, gene=cf$gene_name, distance=cf$distance)
markergenes <- rownames(markers)
markergenes <- rownames(markers)

input_cds <- preprocess_cds(input_cds, num_dim = 54)
input_cds <- align_cds(input_cds, alignment_group = "orig.ident")
input_cds <- reduce_dimension(input_cds, reduction_method = 'UMAP', 
                              preprocess_method = "Aligned")
plot_cells(input_cds, color_cells_by = "orig.ident")
plot_cells(input_cds, color_cells_by = "celltype")
plot_cells(input_cds, color_cells_by = "predicted.id")

input_cds <- cluster_cells(input_cds)
input_cds <- learn_graph(input_cds,close_loop = F)
input_cds <- order_cells(input_cds) 
plot_cells(input_cds, color_cells_by = "orig.ident",cell_size = 1,label_cell_groups=F)
plot_cells(input_cds, color_cells_by = "celltype",cell_size = 1,label_cell_groups=F)
plot_cells(input_cds, color_cells_by = "pseudotime",cell_size = 1,trajectory_graph_color="yellow")
plot_cells(input_cds, color_cells_by = "prediction.score.PT",cell_size = 1,label_cell_groups=T)

input_cds_lin <- input_cds[,is.finite(pseudotime(input_cds))]
plot_accessibility_in_pseudotime(input_cds_lin[c("chr1:100719326-100724136","chr3:171283312-171284976"),],breaks = 8) #VCAM1/TNIK

p <- pseudotime(input_cds)
p <- p[order(p)]
p <- as.data.frame(p)
colnames(p) <- "pseudotime"
#p$pseudotime_order <- c(rep(1:100, each=10)) #8264 cells
#p$pseudotime_order  <- cut(p$pseudotime, breaks = seq(-0.0001,8.303220,length = 11),label=1:10)
p$pseudotime_order <- 1:1000
atacPT <- AddMetaData(atacPT,p)
Idents(atacPT) <- "pseudotime_order"
atacPT@active.ident <- ordered(atacPT@active.ident, levels=1:1000)
atacPT <- ScaleData(atacPT, assay = "chromvar")
#peaks_average <- AverageExpression(atacPT,assays = "chromvar")
peaks_average <- AverageExpression(atacPT,assays = "chromvar",features = c("REL","RELA","FOS::JUN","HNF1A","HNF1B","PPARA::RXRA","NR1H2::RXRA"))
pheatmap::pheatmap(peaks_average[["chromvar"]],cluster_cols = F,cluster_rows = F,scale = "column",border_color = "NA") 
DimPlot(atacPT)

FeaturePlot(atacPT,features = "HNF1A")
