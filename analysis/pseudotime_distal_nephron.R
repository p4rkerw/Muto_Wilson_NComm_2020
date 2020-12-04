# this script will perform pseudotemporal ordering of the distal nephron

# run interactively in local docker container
# navigate to http://localhost:8787
# docker run \
# -v $HOME:$HOME \
# -v /g/diabneph:/home/rstudio/ \
# --rm -p 8787:8787 -e PASSWORD=password p4rkerw/seurat:3.1


library(Seurat)
library(Signac) 
library(monocle3)
library(cicero)
library(BuenColors)
library(EnsDb.Hsapiens.v86)
library(here)
library(dplyr)
library(stringr)
library(data.table)
library(ggplot2)
set.seed(1234)

#Use the object with all version of jasper motifs
atac <- readRDS("cellranger_atac_prep/chromVar.atacAggr_sub97_control.rds")
remotes::install_github('satijalab/seurat-wrappers')
library(SeuratWrappers)


peaks <- data.frame(atac@assays$peaks@counts@Dimnames[1])
colnames(peaks) <- "coord"
peaks.bed <- data.frame(str_split(peaks$coord, pattern=":|-", simplify=TRUE))
colnames(peaks.bed) <- c("chr","start","end")
write.table(peaks.bed, file = "cellranger_atac_prep/peaks.bed", row.names = FALSE, quote = FALSE, col.names = FALSE)

# create cell data set object with cicero constructor
input_cds <- as.cell_data_set(atac)
input_cds <- detect_genes(input_cds)
input_cds <- estimate_size_factors(input_cds)
input_cds <- preprocess_cds(input_cds, num_dim = 50)
input_cds <- align_cds(input_cds, alignment_group = "orig.ident")
input_cds <- reduce_dimension(input_cds, reduction_method = 'UMAP', 
                              preprocess_method = "Aligned")
input_cds <- cluster_cells(input_cds)
plot_cells(input_cds)
input_cds <- learn_graph(input_cds)
# input_cds <- order_cells(input_cds) #Choose the three most distant point from PT_KIM1 in PCT and PST for "start points" of the trajectories. 

# subset for the distal nephron
cds_subset <- input_cds[,colData(input_cds)$celltype %in% c("DCT","CNT","PC")] #Subsetting PT and PT_KIM1 clusters 
cds_subset <- choose_cells(cds_subset) #Select the main clucster and exclude the scattered, low-quality cells.
cds_subset <- cluster_cells(cds_subset)
cds_subset <- learn_graph(cds_subset)
cds_subset <- order_cells(cds_subset) #Choose the most distant point of PC 
plot_cells(cds_subset, color_cells_by = "celltype") +
  ggtitle("Pseudotemporal orderings of distal nephron snATAC dataset")

# create a gene activity matrix for cds_subset
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

cds_subset <- annotate_cds_by_site(cds_subset, gene_annotation_sub)

conns <- fread("analysis_control/ccans/monocle3/ciceroConns.control.allcells.csv") %>%
  dplyr::filter(!is.na(CCAN)) %>% # remove conns that do not belong to a CCAN
  dplyr::select(Peak1, Peak2, coaccess) %>%
  as.data.frame()

unnorm_ga <- build_gene_activity_matrix(cds_subset, conns)

# remove any rows/columns with all zeros
unnorm_ga <- unnorm_ga[!Matrix::rowSums(unnorm_ga) == 0,
                       !Matrix::colSums(unnorm_ga) == 0]

num_genes <- pData(cds_subset)$num_genes_expressed
names(num_genes) <- row.names(pData(cds_subset))

# create a df of activities with colnames as barcodes and rownames as genes
cicero_gene_activities <- normalize_gene_activities(unnorm_ga, num_genes)


# cds_subset@principal_graph_aux@listData$UMAP$root_pr_nodes <- c("Y_32", "Y_65", "Y_66")

fig <- plot_cells(cds_subset,
                    color_cells_by = "celltype",
                    label_groups_by_cluster=FALSE,
                    label_leaves=FALSE,
                    label_branch_points=FALSE,
                    label_roots=FALSE,
                    label_cell_groups=F,show_trajectory_graph=F,group_label_size=3.6) #png 560x430


fig <- plot_cells(cds_subset,
                      color_cells_by = "pseudotime",
                      label_groups_by_cluster=FALSE,
                      label_leaves=FALSE,
                      label_branch_points=FALSE,
                      label_roots=FALSE) #png 560x430

# group cells into 10 bins based on pseudotime
# bins are called "cell_subtype"
cds_subset_lin <- cds_subset[,is.finite(pseudotime(cds_subset))]
pData(cds_subset_lin)$Pseudotime <- pseudotime(cds_subset_lin)
pData(cds_subset_lin)$cell_subtype <- cut(pseudotime(cds_subset_lin), 10)
binned_cds_subset_lin <- aggregate_by_cell_bin(cds_subset_lin, "cell_subtype")

# fit a model to find regions that differ relative to pseudotime
set.seed(1000)
acc_fits <- fit_models(binned_cds_subset_lin[sample(1:nrow(fData(binned_cds_subset_lin)), 1000),], 
                       model_formula_str = "~Pseudotime + num_genes_expressed" )
fit_coefs <- coefficient_table(acc_fits)

# Subset out the differentially accessible sites with respect to Pseudotime
pseudotime_terms <- subset(fit_coefs, term == "Pseudotime" & p_value < .05)
# pseudotime_terms <- readRDS("pseudotime_terms.rds")
head(pseudotime_terms)
df <- as.data.frame(pseudotime_terms)

# take peaks that are differential relative to pseudotime and perform motif enrichment
sub_atac <- subset(atac, idents = c("DCT","CNT","PC"))

# pseudo peaks
pseudo.peak <- df$site_name
# pseudo.peak <- read.csv("pseudo.peak.csv")

# set background peaks
peak.metadata <- GetAssayData(sub_atac, assay = 'peaks', slot = 'meta.features')
background.use <- MatchRegionStats(meta.feature = peak.metadata, regions = pseudo.peak,
                                   features.match = 'GC.percent', n = 10000)

# find enriched motifs
enriched.motifs <- FindMotifs(
  object = sub_atac,
  features = pseudo.peak,
  background = background.use
)

# filter motifs
enriched.motifs.filter <- dplyr::filter(enriched.motifs, pvalue < 0.05) %>%
  dplyr::arrange(-fold.enrichment)

# annotate the motifs with gene names
# load the ensembl annotations
gene.coords <- genes(EnsDb.Hsapiens.v86, filter = ~ gene_biotype == "protein_coding", columns = "gene_name")
seqlevelsStyle(gene.coords) <- 'UCSC'
gene.coords <- keepStandardChromosomes(gene.coords, pruning.mode = 'coarse')
gene.coords <- Extend(gene.coords, upstream=50000, downstream=50000)

pseudo_peaks.bed <- str_split(pseudo.peak, pattern = ":|-", simplify = TRUE) 
colnames(pseudo_peaks.bed) <- c("chrom","start","end")
pseudo_peaks.gr <- makeGRangesFromDataFrame(data.frame(pseudo_peaks.bed))
BiocManager::install("plyranges")
library(plyranges)

# annotate the pseudotemporal diff accessible regions with overlapping ensembl genes
anno_atac <- join_overlap_left(pseudo_peaks.gr, gene.coords)
anno_atac <- sortSeqlevels(anno_atac)
anno_atac <- sort(anno_atac)


# perform the same analysis in the corresponding snRNA celltypes
rna <- readRDS("cellranger_rna_prep/rnaAggr_control.rds")

# create cell data set object with cicero constructor
rna_input_cds <- as.cell_data_set(rna)
rna_input_cds <- detect_genes(rna_input_cds)
rna_input_cds <- estimate_size_factors(rna_input_cds)
rna_input_cds <- preprocess_cds(rna_input_cds, num_dim = 50)
rna_input_cds <- align_cds(rna_input_cds, alignment_group = "orig.ident")
rna_input_cds <- reduce_dimension(rna_input_cds, reduction_method = 'UMAP', 
                              preprocess_method = "Aligned")
rna_input_cds <- cluster_cells(rna_input_cds)
plot_cells(rna_input_cds, color_cells_by = "celltype")

rna_cds_subset <- rna_input_cds[,colData(rna_input_cds)$celltype %in% c("DCT1","DCT2","CNT","PC")] #Subsetting PT and PT_KIM1 clusters 
rna_cds_subset <- choose_cells(rna_cds_subset) #Select the main clucster and exclude the scattered, low-quality cells.
rna_cds_subset <- cluster_cells(rna_cds_subset)
plot_cells(rna_cds_subset, color_cells_by = "celltype")
rna_cds_subset <- learn_graph(rna_cds_subset)
rna_cds_subset <- order_cells(rna_cds_subset) #Choose the most distant point of DCT

fig <- plot_cells(rna_cds_subset,
                  color_cells_by = "celltype",
                  label_groups_by_cluster=FALSE,
                  label_leaves=FALSE,
                  label_branch_points=FALSE,
                  label_roots=FALSE) +
ggtitle("Pseudotemporal ordering of distal nephron snRNA dataset")

fig <- plot_cells(rna_cds_subset,
                  color_cells_by = "pseudotime",
                  label_groups_by_cluster=FALSE,
                  label_leaves=FALSE,
                  label_branch_points=FALSE,
                  label_roots=FALSE) #png 560x430

# find genes that differ rel to pseudotime
cds_pr_test_res <- graph_test(rna_cds_subset, neighbor_graph="principal_graph", cores=4)
pr_deg_ids <- row.names(subset(cds_pr_test_res, q_value < 0.01))
gene_module_df <- find_gene_modules(rna_cds_subset[pr_deg_ids,], resolution=0.001)

cell_group_df <- tibble::tibble(cell=row.names(colData(rna_cds_subset)), 
                                cell_group=colData(rna_cds_subset)$celltype)
agg_mat <- aggregate_gene_expression(rna_cds_subset, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
pheatmap::pheatmap(agg_mat,
                   scale="column", clustering_method="ward.D2")

rowData(rna_cds_subset)$gene_short_name <- rownames(rna_cds_subset) # there is a bug in seurat import fn
plot_cells(rna_cds_subset,
           genes=gene_module_df %>% filter(module %in% c(6,8,12)),
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)

# intersect 
module6_genes <- gene_module_df %>% filter(module %in% 6)
int <- intersect(module6_genes$id, anno_atac$gene_name)
int2 <- intersect(enriched.motifs.filter$motif.name, module6_genes$id)


