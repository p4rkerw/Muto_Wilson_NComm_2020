library(Signac) # 0.2.1
library(Seurat) # 3.0.2
library(cicero) # 1.3.4
library(monocle3) # 3.0.2
library(dplyr)
library(here)
library(openxlsx)
library(rtracklayer)
library(data.table)
library(tibble)
library(ggplot2)

# input Seurat object
PrepareInputCDS <- function(seurat_obj) {
  #seurat_obj <- subset(seurat_agg, ident=celltype)
  DefaultAssay(seurat_obj) <- "peaks"
  count_data <- GetAssayData(seurat_obj, slots = "counts")
  summ <- summary(count_data)
  summ_frame <- data.frame(peak = rownames(count_data)[summ$i],
                           cell.id = colnames(count_data)[summ$j],
                           count = summ$x)
  summ_frame <- RenamePeaks(summ_frame) # convert peaks to cicero naming convention
  
  # create cell data set object with cicero constructor
  input_cds <- make_atac_cds(summ_frame, binarize = TRUE)
  
  # input_cds <- input_cds[Matrix::rowSums(exprs(input_cds)) != 0,] 
  set.seed(2017)
  input_cds <- detect_genes(input_cds)
  input_cds <- estimate_size_factors(input_cds)
  input_cds <- preprocess_cds(input_cds, method="LSI")
  input_cds <- reduce_dimension(input_cds, reduction_method="UMAP", preprocess_method="LSI")
  umap_coords <- reducedDims(input_cds)$UMAP
  return(input_cds)
}

RenamePeaks <- function(summ_frame) {
  require(stringr)
  
  peaks <- summ_frame$peak
  peaks <- str_split(peaks, pattern = ":", simplify = TRUE)
  peaks <- cbind(peaks[, 1], str_split(peaks[, 2], pattern = "-", simplify = TRUE))
  peaks_formatted <- paste(peaks[, 1], peaks[, 2], peaks[, 3], sep = "_")
  summ_frame$peak <- peaks_formatted
  return(summ_frame)
}

# load the aggregated atac object
atacAggr <- readRDS("cellranger_atac_prep/chromVar.atacAggr_sub97_control.rds")

# load gene annotation and format
temp <- tempfile()
download.file("ftp://ftp.ensembl.org/pub/release-99/gtf/homo_sapiens/Homo_sapiens.GRCh38.99.gtf.gz", temp)
gene_anno <- rtracklayer::readGFF(temp)
unlink(temp)
gene_anno$chromosome <- paste0("chr", gene_anno$seqid)
gene_anno$gene <- gene_anno$gene_id
gene_anno$transcript <- gene_anno$transcript_id
gene_anno$symbol <- gene_anno$gene_name

#### Add a column for the pData table indicating the gene if a peak is a promoter ####
# Create a gene annotation set that only marks the transcription start sites of 
# the genes. We use this as a proxy for promoters.
# To do this we need the first exon of each transcript
pos <- subset(gene_anno, strand == "+")
pos <- pos[order(pos$start),] 
# remove all but the first exons per transcript
pos <- pos[!duplicated(pos$transcript),] 
# make a 1 base pair marker of the TSS
pos$end <- pos$start + 1 

neg <- subset(gene_anno, strand == "-")
neg <- neg[order(neg$start, decreasing = TRUE),] 
# remove all but the first exons per transcript
neg <- neg[!duplicated(neg$transcript),] 
neg$start <- neg$end - 1
gene_annotation_sub <- rbind(pos, neg)

# Make a subset of the TSS annotation columns containing just the coordinates 
# and the gene name
gene_annotation_sub <- gene_annotation_sub[,c("chromosome", "start", "end", "symbol")]

# Rename the gene symbol column to "gene"
names(gene_annotation_sub)[4] <- "gene"

input_cds <- PrepareInputCDS(atacAggr)

input_cds <- annotate_cds_by_site(input_cds, gene_annotation_sub)

conns <- fread("analysis_control/ccans/monocle3/ciceroConns.control.allcells.csv") %>%
  dplyr::filter(!is.na(CCAN)) %>% # remove conns that do not belong to a CCAN
  dplyr::select(Peak1, Peak2, coaccess) %>%
  as.data.frame()

#### Generate gene activity scores ####
# generate unnormalized gene activity matrix
unnorm_ga <- build_gene_activity_matrix(input_cds, conns)

# remove any rows/columns with all zeros
unnorm_ga <- unnorm_ga[!Matrix::rowSums(unnorm_ga) == 0,
                       !Matrix::colSums(unnorm_ga) == 0]

num_genes <- pData(input_cds)$num_genes_expressed
names(num_genes) <- row.names(pData(input_cds))

# create a df of activities with colnames as barcodes and rownames as genes
cicero_gene_activities <- normalize_gene_activities(unnorm_ga, num_genes)

# create a cicero gene activity matrix for each cell type
Idents(atacAggr) <- "highres.predicted.id"
atacAggr <- RenameIdents(atacAggr, "PT_KIM1" = "PT_VCAM1")
Idents(atacAggr) -> atacAggr@meta.data$highres.predicted.id
celltype <- levels(atacAggr)
cell_cicero_gene_act.ls <- lapply(celltype, function(ident){
  meta <- data.frame(barcode=rownames(atacAggr@meta.data), 
                     celltype=atacAggr@meta.data$highres.predicted.id)
  cell_bc <- dplyr::filter(meta, celltype == ident)
  cicero_gene_act <- cicero_gene_activities[,cell_bc$barcode]
  mean_cicero_act <- rowMeans(cicero_gene_act) %>% as.data.frame()
  df <- data.frame(gene = rownames(mean_cicero_act), mean_cicero_act=mean_cicero_act[,1])
  return(df)
}) 
names(cell_cicero_gene_act.ls) <- celltype

# find all cell-specific markers
rnaAggr <- readRDS("cellranger_rna_prep/rnaAggr_control.rds")
rnaAggr <- NormalizeData(rnaAggr, assay = 'RNA')
cell_deg.ls <- lapply(celltype, function(ident){
  df <- read.xlsx("analysis_control/deg.celltype.control.xlsx", sheet = ident, rowNames = TRUE)
  df$celltype <- ident
  df <- rownames_to_column(df, var = "gene")
  return(df)
})
names(cell_deg.ls) <- celltype

# join the cicero activity and gene expression dataframe
merge.ls <- lapply(celltype, function(ident){
  cicero_act <- cell_cicero_gene_act.ls[[ident]]
  gene_exp <- cell_deg.ls[[ident]]
  mat <- tryCatch(full_join(cicero_act, gene_exp, by = "gene"),
                  error=function(e) NULL)  
  # remove rows that do not have gene expression or cicero activity
  mat <- dplyr::filter(mat, !is.na(avg_logFC), !is.na(mean_cicero_act))
  return(mat)
})

# aggregate the merged df
df <- bind_rows(merge.ls)
df <- dplyr::arrange(df, gene, celltype)
df$celltype <- factor(df$celltype, levels=celltype)
write.csv(df, file = "analysis_control/corr_cicero_gene_exp.csv")

# filter the motif-gene df for pairs with significant changes in gene expression and chromvar activity
df.f <- dplyr::filter(df, p_val_adj < 0.05) 
df.f <- dplyr::mutate(df.f, combo = paste0(gene,"_",celltype))
df.f <- dplyr::distinct(df.f, combo, .keep_all=TRUE)


# calculate pearson r2 for all cicero-gene combos
pearson <- cor.test(df.f$mean_cicero_act, df.f$avg_logFC, method="pearson", conf.level=0.95)
max_cicero <- max(abs(df.f$mean_cicero_act)) * 1.1
max_exp <- max(abs(df.f$avg_logFC)) * 1.1

# draw a global plot of cicero activity vs. expression
p1 <- ggplot(df.f, aes(x=mean_cicero_act, y=avg_logFC)) +
  geom_smooth(method="lm", color = "darkgray") +
  geom_point(aes(color=celltype)) +
  xlab("Cicero activity") +
  ylab("Gene expression (avg_logFC)") +
  ggtitle("Cicero Gene Activity vs. Gene Expression for all celltypes",
          subtitle=paste0("Pearson r^2=",
                          signif(pearson$estimate,2),
                          " pval=",signif(pearson$p.val,2))) +
  xlim(c(0,max_cicero)) +
  ylim(c(-max_exp,max_exp)) +
  theme_minimal() +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) 
p1

# calculate a pearson correlation between cicero activity and gene expression for every detected gene
library(ggpmisc)
df.pearson <- lapply(df.f$gene, function(gene_sel) {
  df <- dplyr::filter(df.f, gene == gene_sel)
  pearson <- tryCatch(cor.test(df$mean_cicero_act, df$avg_logFC,
                               method="pearson", conf.level=0.95),
                      error=function(e) NULL)
  df$corr <- tryCatch(signif(pearson$estimate,2), error=function(e) NULL)
  df$pval <- tryCatch(signif(pearson$p.value,2), error=function(e) NULL)
  df$max_exp <- max_exp
  df$max_chrom <-max_chrom
  df$num_celltypes <- length(unique(df$celltype))
  return(df)
})
df.final <- bind_rows(df.pearson)
df.final <- dplyr::distinct(df.final, combo, .keep_all=TRUE)

# filter for significant correlations
df.sig <- dplyr::filter(df.final, pval < 0.05) %>%
  arrange(-num_celltypes, corr)
library(ggrepel)
gene_sel = "NR3C2"
gene.f <- dplyr::filter(df.sig, gene==gene_sel)
p1 <- ggplot(gene.f, aes(x=mean_cicero_act, y=avg_logFC)) +
  geom_smooth(method="lm", color = "darkgray") +
  geom_point(aes(color=celltype)) +
  geom_text_repel(aes(label=celltype)) +
  xlab("Cicero activity") +
  ylab("Gene expression (avg_logFC)") +
  ggtitle(gene_sel,
          subtitle=paste0("Pearson r^2=",
                          unique(df.sig$corr),
                          " pval=",unique(df.sig$pval))) +
  theme_minimal() +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) 
p1

