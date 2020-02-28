# this script will generate cis-coaccessibility networks for individual celltypes. 
library(Signac) # 0.1.5
library(Seurat) # 3.0.2
library(cicero) # 1.2
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
  input_cds <- detectGenes(input_cds)
  print(colnames(pData(input_cds)))
  input_cds <- estimateSizeFactors(input_cds)
  
  # *** if you are using Monocle 3, you need to run the following line as well!
  #input_cds <- preprocessCDS(input_cds, norm_method = "none")
  input_cds <- reduceDimension(input_cds, max_components = 2, num_dim=6,
                               reduction_method = 'tSNE', norm_method = "none")
  
  tsne_coords <- t(reducedDimA(input_cds)) # dimensional reduction coords could be extracted from seurat object
  row.names(tsne_coords) <- row.names(pData(input_cds))
  cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = tsne_coords) 
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

Get_Ccans <- function(clusterID, seurat_agg) {
  print(paste0("Calculating CCAN for: ",clusterID))
  atacAggr <- seurat_agg
  seuratSub <- subset(atacAggr, ident = clusterID) # create a subset
  
  # convert seurat objects into cicero cell datasets in preparation for detecting cicero connections
  ciceroCds <- PrepareCiceroCDS(seuratSub)
  
  # generate disease-specific CCANS for all chromsomes of a particular celltype
  # FindCiceroConns can also take specific seqLevels eg.
  # seqLevels <- paste0("chr",c(21,22))
  # FindCiceroConns(cds, seqLevels)
  conns <- FindCiceroConns(ciceroCds)
  return(conns)
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
  write.csv(list.ccan[[x]], file = paste0("analysis_control/ccans/ciceroConns.control.",x,".csv", row.names = TRUE))
})




