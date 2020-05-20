# this script will generate cis-coaccessibility networks for individual celltypes. 
library(Signac) # 0.2.1
library(Seurat) # 3.0.2
library(cicero) # 1.3.4
library(monocle3) # 3.0.2
library(dplyr)
library(here)
library(openxlsx)

# input Seurat object
PrepareCiceroCDS <- function(seurat_obj) {
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
  cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates=umap_coords)
  
  return(cicero_cds)
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

FindCiceroConns <- function(cds, chrom = NULL) {
  # to generate the grch38 contigs execute the following code in linux
  # python getContigLengths.py /home/parkerw/reference/refdata-cellranger-GRCh38-3.0.0/fasta/genome.fa temp.txt
  # cat temp.txt | sed 's/>//g' | sed 's/GRCh38_//g' > contigLengths.txt
  # rm temp.txt
  
  # Note make sure to use the correct reference 
  # hg19 coordinates are available preloaded in the package but GRCh38 has to be created (see above)
  # data("human.hg19.genome")
  contigs <- read.table("github_repository/healthyKidney/utility_scripts/contigLengths.txt")
  contigs$V1 <- paste0("chr",contigs$V1)
  contigs <- select(contigs, c("V1","V5"))
  
  # select chromosomes to run cicero on
  levels <- paste0("chr",c(seq(1,22),"X","Y"))
  contigs <- subset(contigs, V1 %in% levels)
  
  # if a specific chromosome is specified, subset the contigs and only run that chromosome
  if(!is.null(chrom)) {
    contigs <- subset(contigs, V1 %in% chrom)
  }
  
  # can the contig region be limited to 1Mb up- and downstream of gene of interest to increase speed?
  # create a subset by input chromosome
  conns <- run_cicero(cds, contigs) 
  return(conns)
}

Get_Ccans <- function(clusterID=NULL, seurat_agg) {
  if(!is.null(clusterID)) {
    print(paste0("Subsetting seurat object for: ",clusterID))
    seurat_agg <- subset(seurat_agg, ident = clusterID) # create a subset
  } 
  # convert seurat objects into cicero cell datasets in preparation for detecting cicero connections
  print("Preparing Cicero CDS")
  ciceroCds <- PrepareCiceroCDS(seurat_agg)
  
  # generate disease-specific CCANS for all chromsomes of a particular celltype
  # FindCiceroConns can also take specific seqLevels eg.
  # seqLevels <- paste0("chr",c(21,22))
  # FindCiceroConns(cds, seqLevels)
  print("Finding Cicero connections")
  conns <- FindCiceroConns(ciceroCds)
  
  CCAN_assigns <- generate_ccans(conns)

  # create a column that identifies which connections belong to a CCAN
  ccan1 <- left_join(conns, CCAN_assigns, by=c("Peak1" = "Peak"), all.x=TRUE)
  colnames(ccan1)[4] <- "CCAN1"
  ccan2 <- left_join(conns, CCAN_assigns, by=c("Peak2" = "Peak"), all.x=TRUE)
  colnames(ccan2)[4] <- "CCAN2"
  df <- cbind(ccan1, CCAN2=ccan2$CCAN2) %>%
    dplyr::mutate(CCAN = ifelse(CCAN1 == CCAN2, CCAN1, 0)) %>%
    dplyr::select(-CCAN1, -CCAN2)
  
  return(df)
}


####################################################################
# load aggregated snATACseq dataset and select cluster of interest
atacAggr <- readRDS("cellranger_atac_prep/atacAggr_sub97_control.rds")

# subset the snATACseq object by disease state within celltype of interest and find ccans
Idents(atacAggr) <- "celltype"
list.ccan <- lapply(levels(atacAggr), function(ident) {Get_Ccans(ident, seurat_agg = atacAggr)})

# write to file
dir.create("analysis_control/ccans", showWarnings = FALSE)
names(list.ccan) <- levels(atacAggr)
lapply(names(list.ccan), function(x) {
  fwrite(list.ccan[[x]], file = paste0("analysis_control/ccans/ciceroConns.control.",x,".csv"), row.names = TRUE)
})

# calculate a global CCAN for all celltypes
ccan <- Get_Ccans(seurat_agg = atacAggr)
fwrite(ccan, file = "analysis_control/ccans/ciceroConns.control.allcells.csv", row.names = TRUE)

# identify which peaks belong to a CCAN
CCAN_assigns <- generate_ccans(ccan)

# create a column that identifies which connections belong to a CCAN
ccan1 <- left_join(ccan, CCAN_assigns, by=c("Peak1" = "Peak"), all.x=TRUE)
colnames(ccan1)[4] <- "CCAN1"
ccan2 <- left_join(ccan, CCAN_assigns, by=c("Peak2" = "Peak"), all.x=TRUE)
colnames(ccan2)[4] <- "CCAN2"
df <- cbind(ccan1, CCAN2=ccan2$CCAN2) %>%
  dplyr::mutate(CCAN = ifelse(CCAN1 == CCAN2, CCAN1, 0)) %>%
  dplyr::select(-CCAN1, -CCAN2)

fwrite(df, file = "analysis_control/ccans/ciceroConns.control.allcells.csv", row.names = TRUE)
