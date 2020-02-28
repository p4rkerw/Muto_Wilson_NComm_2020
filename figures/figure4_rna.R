library(Seurat)
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

genes <- c("VCAM1","TPM1","SLC5A12","SLC4A4")
lineage_cds <- cds_subset[rowData(cds_subset)$gene_short_name %in% genes,]
fig4d <- plot_genes_in_pseudotime(lineage_cds,
                         color_cells_by="celltype",
                         min_expr=1,
                         cell_size=0.1,
                         trend_formula = "~ splines::ns(pseudotime, df=10)",
                         panel_order = c("VCAM1","TPM1","SLC5A12","SLC4A4")
                         ) #png 500x670

genes <- c("VCAM1","TPM1")
lineage_cds <- cds_subset[rowData(cds_subset)$gene_short_name %in% genes,]
fig4d_1 <- plot_genes_in_pseudotime(lineage_cds,
                                  color_cells_by="celltype",
                                  min_expr=1,
                                  cell_size=0.1,
                                  trend_formula = "~ splines::ns(pseudotime, df=10)",
                                  panel_order = c("VCAM1","TPM1")
) #png 500x670

genes <- c("SLC5A12","SLC4A4")
lineage_cds <- cds_subset[rowData(cds_subset)$gene_short_name %in% genes,]
fig4d_2 <- plot_genes_in_pseudotime(lineage_cds,
                                  color_cells_by="celltype",
                                  min_expr=1,
                                  cell_size=0.1,
                                  trend_formula = "~ splines::ns(pseudotime, df=10)",
                                  panel_order = c("SLC5A12","SLC4A4")
) #png 500x670