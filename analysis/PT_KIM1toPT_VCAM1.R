#snRNA-seq dataset
Idents(rnaAggr) <- "seurat_clusters"
new.cluster.ids <- c("PT","DCT1","TAL","CNT","PT",
                     "TAL","ICA","TAL","PC","ENDO",
                     "PEC","DCT2","PODO","PT_VCAM1","ICB",
                     "ENDO","MES","FIB","ENDO","LEUK")
names(new.cluster.ids) <- levels(rnaAggr)
rnaAggr <- RenameIdents(rnaAggr, new.cluster.ids)


# reorder the idents and save celltype annotations in celltype slot metadata
levels(rnaAggr) <- c("PT","PT_VCAM1","PEC","TAL","DCT1",
                     "DCT2","CNT","PC","ICA","ICB",
                     "PODO","ENDO","MES","FIB","LEUK")
rnaAggr@meta.data$celltype <- rnaAggr@active.ident
DimPlot(rnaAggr,label = T,reduction = "umap")
#saveRDS(rnaAggr,"cellranger_rna_prep/rnaAggr_control.rds")

#snATAC-seq dataset
Idents(sub_atac) <- "seurat_clusters"
new.cluster.ids <- c("PCT","TAL","PST","PCT","TAL",
                     "DCT","TAL","PST","PC","CNT",
                     "ENDO","PT_VCAM1","ICB","ICA","TAL",
                     "PEC","DCT","MES_FIB","LEUK","PODO")

names(new.cluster.ids) <- levels(sub_atac)
sub_atac <- RenameIdents(sub_atac, new.cluster.ids)
levels(sub_atac) <- c("PCT","PST","PT_VCAM1","PEC","TAL",
                      "DCT","CNT","PC","ICA","ICB",
                      "PODO","ENDO","MES_FIB","LEUK")
sub_atac@meta.data$celltype <- sub_atac@active.ident
DimPlot(sub_atac,label = T,reduction = "umap")
#saveRDS(sub_atac,"cellranger_atac_prep/atacAggr_sub97_control.rds")

