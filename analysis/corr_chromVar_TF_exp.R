# this script will correlate chromVAR transcription factor activity with TF expression

library(Seurat)
library(Signac) # 0.2.5
library(Seurat) # 3.1.5
library(EnsDb.Hsapiens.v86) 
library(GenomicRanges) # 1.38.0
library(harmony) # 1.0
library(BSgenome.Hsapiens.UCSC.hg38) # 1.4.1
library(future) # 1.17.0
library(ggplot2)
library(dplyr) # 1.0.0
library(here) # 0.1
library(stringr)

atacAggr <- readRDS("cellranger_atac_prep/chromVar.atacAggr_sub97_control.rds")

rnaAggr <- readRDS("cellranger_rna_prep/rnaAggr_control.rds")
rnaAggr <- NormalizeData(rnaAggr, assay = 'RNA')

transfer.anchors <- FindTransferAnchors(
  reference = rnaAggr,
  query = atacAggr,
  reference.assay = 'RNA',
  query.assay = 'RNA',
  reduction = 'cca'
)

refdata <- GetAssayData(
  object = rnaAggr, 
  assay = "RNA", 
  slot = "data"
)

imputation <- TransferData(
  anchorset = transfer.anchors, 
  refdata = refdata, 
  weight.reduction = atacAggr[["pca"]] 
)

atacAggr[["IMPRNA"]] <- imputation

saveRDS(transfer.anchors, file = here("cellranger_atac_prep/transfer_anchors_control.rds"))
rm(transfer.anchors)

# saveRDS(atacAggr, file = here("cellranger_atac_prep/imputed.chromVar.atacAggr_sub97_control.rds"))
atacAggr <- readRDS("cellranger_atac_prep/imputed.chromVar.atacAggr_sub97_control.rds")
Idents(atacAggr) <- "highres.predicted.id"
atacAggr <- RenameIdents(atacAggr, "PT_KIM1" = "PT_VCAM1")

# grab celltype annotation
celltype <- atacAggr@meta.data$predicted.id

# create dotplot of chromVar activity vs. imputed gene expression
chromvar <- GetAssayData(atacAggr, assay="chromvar", slot="data")
motifs <- atacAggr$peaks@misc$motif@motif.names

# create a data frame of motifs and corresponding gene names
motif_names <- data.frame(genes=unlist(motifs)) %>%
  rownames_to_column(var="motif") %>%
  arrange(genes)

# format the motif_names df so each motif corresponds to a single gene
# the JASPAR database has multiple genes associated with each motif separated by "::"
motif_names <- dplyr::mutate(motif_names, gene = str_split(genes, pattern="::")) %>%
  unnest(gene) %>%
  dplyr::select(motif, gene) %>%
  distinct() %>%
  arrange(gene) %>%
  dplyr::filter(!str_detect(gene, pattern="var.")) %>% # remove any genes associated with multiple variants
  as.data.frame()

# fix the one occurrence of EWSR1-FLI1
motif_names$gene <- str_replace(motif_names$gene, pattern="EWSR1-FLI1", replacement="EWSR1")

# use the findmarkers function to calculate gene expression and activity for every motif-gene combination
motif.ls = motif_names$motif

df.ls <- lapply(motif.ls, function(motif) {
    print(motif)
    gene <- motif_names[motif_names$motif == motif,]$gene
    print(gene)
    DefaultAssay(atacAggr) <- "chromvar"
    aver_chromvar <- FindAllMarkers(atacAggr, min.pct=0, logfc.threshold=0, features=motif,
                                    verbose=FALSE)
    
    # compute the average expression from the rna object
    aver_exp <- FindAllMarkers(rnaAggr, assays="RNA",min.pct=0, logfc.threshold=0, features=gene, verbose=FALSE)
    
    # join the matrices and return null if gene expression is not detected for all cell types
    mat <- tryCatch(full_join(aver_chromvar, aver_exp, by = "cluster"),
                              error=function(e) NULL)
    if(!is.null(mat)) {
      df <- data.frame(chromvar=mat$avg_logFC.x, rna=mat$avg_logFC.y, celltype=mat$cluster,
                     motif=mat$gene.x, gene=mat$gene.y, chrom_pval=mat$p_val_adj.x, gene_pval=mat$p_val_adj.y)
    } else {
      return(NULL)
    }
})

# aggregate the motif-gene df
df <- bind_rows(df.ls)
df <- dplyr::arrange(df, gene, motif)
write.csv(df, file = "analysis_control/corr_chromVar_TF_exp.csv")

# filter the motif-gene df for pairs with significant changes in gene expression and chromvar activity
df.f <- dplyr::filter(df, chrom_pval < 0.05, gene_pval < 0.05)
df.f <- na.omit(df.f)

# create a column for each unique gene-motif combo
df.f <- dplyr::mutate(df.f, gene_motif = paste0(gene,"_",motif))

# compute a pearson r2 and pval for each motif-gene combo across all celltypes
combos <- unique(df.f$gene_motif)
library(ggpmisc)
df.pearson <- lapply(combos, function(combo) {
  require(ggrepel)

  df <- dplyr::filter(df.f, gene_motif == combo)
  max_exp <- max(abs(df$rna), na.rm=TRUE) + 1
  max_chrom <- max(abs(df$chromvar), na.rm=TRUE) + 1
  pearson <- tryCatch(cor.test(df$chromvar, df$rna, method="pearson", conf.level=0.95), error=function(e) NULL)
  df$corr <- tryCatch(signif(pearson$estimate,2), error=function(e) NULL)
  df$pval <- tryCatch(signif(pearson$p.value,2), error=function(e) NULL)
  df$max_exp <- max_exp
  df$max_chrom <-max_chrom
  df$num_celltypes <- length(unique(df$celltype))
  return(df)
})

# aggregate a final df with all the available pearson R2 coefficients and pvals
df.final <- bind_rows(df.pearson)
df.final <- dplyr::distinct(df.final, celltype, gene_motif, .keep_all = TRUE)
df.sig <- dplyr::filter(df.final, pval < 0.05) %>%
  arrange(-num_celltypes, -corr)

# calculate pearson r2 for all tf-gene combos
pearson <- cor.test(df.sig$chromvar, df.sig$rna, method="pearson", conf.level=0.95)
max_chrom <- max(abs(df.sig$max_chrom)) + 1
max_exp <- max(abs(df.sig$max_exp)) + 1

p1 <- ggplot(df.sig, aes(x=chromvar, y=rna)) +
  geom_smooth(method="lm", color = "darkgray") +
  geom_point(aes(color=celltype)) +
  xlab("chromVAR activity (avg_logFC)") +
  ylab("Gene expression (avg_logFC)") +
  ggtitle("Significant gene-motif combos for all celltypes", subtitle=paste0("Pearson r^2=",signif(pearson$estimate,2)," pval=",signif(pearson$p.val,2))) +
  xlim(c(-max_chrom,max_chrom)) +
  ylim(c(-max_exp,max_exp)) +
  theme_minimal() +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) 
p1


# identify TF that are significantly correlated with their chromvar activity
sig_pairs <- unique(df.sig$gene_motif)
total_pairs <- unique(df.final$gene_motif)

# identify TF that are positively (or negatively) correlated with their chromvar activity
df.pos <- dplyr::filter(df.sig, corr > 0)
sig_pos <- unique(df.pos$gene_motif)
pearson <- cor.test(df.pos$chromvar, df.pos$rna, method="pearson", conf.level=0.95)
max_chrom <- max(abs(df.pos$max_chrom)) + 1
max_exp <- max(abs(df.pos$max_exp)) + 1

p1 <- ggplot(df.pos, aes(x=chromvar, y=rna)) +
  geom_smooth(method="lm", color = "darkgray") +
  geom_point(aes(color=celltype)) +
  xlab("chromVAR activity (avg_logFC)") +
  ylab("Gene expression (avg_logFC)") +
  ggtitle("Positive correlation gene-motif combos for all celltypes", 
          subtitle=paste0("Pearson r^2=",signif(pearson$estimate,2)," pval=",signif(pearson$p.val,2))) +
  xlim(c(-max_chrom,max_chrom)) +
  ylim(c(-max_exp,max_exp)) +
  theme_minimal() +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) 
p1



df.neg <- dplyr::filter(df.sig, corr < 0)
sig_neg <- unique(df.neg$gene_motif)
pearson <- cor.test(df.neg$chromvar, df.neg$rna, method="pearson", conf.level=0.95)
max_chrom <- max(abs(df.neg$max_chrom)) + 1
max_exp <- max(abs(df.neg$max_exp)) + 1

p2 <- ggplot(df.neg, aes(x=chromvar, y=rna)) +
  geom_smooth(method="lm", color = "darkgray") +
  geom_point(aes(color=celltype)) +
  xlab("chromVAR activity (avg_logFC)") +
  ylab("Gene expression (avg_logFC)") +
  ggtitle("Negative correlation gene-motif combos for all celltypes", 
          subtitle=paste0("Pearson r^2=",signif(pearson$estimate,2)," pval=",signif(pearson$p.val,2))) +
  xlim(c(-max_chrom,max_chrom)) +
  ylim(c(-max_exp,max_exp)) +
  theme_minimal() +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) 
p2


toplot <- dplyr::filter(df.final, gene_motif == "ZEB1_MA0103.3")
max_chrom <- unique(toplot$max_chrom)
max_exp <- unique(toplot$max_exp)
p1 <- ggplot(toplot, aes(x=chromvar, y=rna)) +
  geom_smooth(method="lm", color = "darkgray") +
  geom_point(aes(color=celltype)) +
  geom_text_repel(aes(label=celltype)) +
  xlab("chromVAR activity (avg_logFC)") +
  ylab("Gene expression (avg_logFC)") +
  ggtitle(unique(toplot$gene), subtitle=paste0("Pearson r^2=",unique(toplot$cor)," pval=",unique(toplot$pval))) +
  xlim(c(-max_chrom,max_chrom)) +
  ylim(c(-max_exp,max_exp)) +
  theme_minimal() +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) 
p1

