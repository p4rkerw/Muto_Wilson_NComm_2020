# this script will generate ucsc formatted bed and interact tracks for visualization
# of DAR and CCAN
library(openxlsx)
library(dplyr)
library(tibble)
library(stringr)
library(GenomicRanges)
library(data.table)

outDir <- "analysis_control/ucsc_tracks/"
dir.create(outDir)

celltypes <- getSheetNames("analysis_control/dar.celltype.control.xlsx")

# read in celltype DAR and create bed tracks to visualize in UCSC
lapply(celltypes, function(celltype) {
  clusterID <- celltype
  dar <- read.xlsx("analysis_control/dar.celltype.control.xlsx", sheet = clusterID, rowNames=TRUE) %>%
    rownames_to_column(var="coord") %>%
    dplyr::mutate(chrom = str_split(coord, pattern=":|-", simplify=TRUE)[,1]) %>%
    dplyr::mutate(chromStart = str_split(coord, pattern=":|-", simplify=TRUE)[,2]) %>%
    dplyr::mutate(chromEnd = str_split(coord, pattern=":|-", simplify=TRUE)[,3]) %>%
    # dplyr::filter(gene == gene.sel) %>%
    dplyr::select(chrom, chromStart, chromEnd, avg_logFC) %>%
    dplyr::mutate(avg_logFC = format(round(avg_logFC, 2), nsmall=2)) # limit to two decimal places
  colnames(dar)[1] <- "#chrom"
  trackName = paste0(clusterID,"_DAR")
  trackDescription = "kidney celltype differentially accessible chromatin"
  comment1 = paste0("track type=bed name=", trackName, " description=\"", trackDescription, "\" color=0,128,0,")
  comment2 = paste0("browser position chr20:44,080,406-44,753,709") 
  fileConn <- file(paste0(outDir, clusterID, ".dar.bed"), "w")
  writeLines(comment1, fileConn)
  writeLines(comment2, fileConn)
  write.table(dar, file=fileConn, sep = "\t", row.names = FALSE, quote = FALSE)
  close(fileConn)
})


# create celltype interaction tracks from cicero CCAN xlsx files
ccanFiles <- list.files("analysis_control/ccans/", pattern = "ciceroConns")
ccanFilePaths <- paste0("analysis_control/ccans/", ccanFiles)
  
lapply(ccanFilePaths, function(filepath) {

  ccan <- fread(filepath)
  clusterID <- str_split(basename(filepath), pattern = "[.]", simplify=TRUE)[,3]
  
  ccan.df <- ccan %>%
    dplyr::filter(abs(coaccess) > 0.2) %>%
    dplyr::mutate(score = as.integer(1000 * abs(coaccess))) %>% # convert cicero coaccess to ucsc score format
    dplyr::mutate(seqnames = str_split(Peak1, pattern="_", simplify=TRUE)[,1]) %>%
    dplyr::mutate(start = str_split(Peak1, pattern="_", simplify=TRUE)[,2]) %>%
    dplyr::mutate(end = str_split(Peak1, pattern="_", simplify=TRUE)[,3]) %>%
    dplyr::mutate(targetChrom = str_split(Peak2, pattern="_", simplify=TRUE)[,1]) %>%
    dplyr::mutate(targetStart = str_split(Peak2, pattern="_", simplify=TRUE)[,2]) %>%
    dplyr::mutate(targetEnd = str_split(Peak2, pattern="_", simplify=TRUE)[,3]) %>%
    dplyr::select(seqnames, start, end, targetChrom, targetStart, targetEnd, score) 
    
  ccan.interact = data.frame(
    chrom=ccan.df[,"seqnames"],
    chromStart=ccan.df[,"start"],
    chromEnd=ccan.df[,"end"],
    name = "CCAN",
    score = ccan.df$score,
    value = ccan.df$score,
    exp = "kidney_chromatin_interactions",
    color = "#0000FF", # make all the arcs the same color
    sourceChrom=ccan.df[,"seqnames"],
    sourceStart=ccan.df[,"start"],
    sourceEnd=ccan.df[,"end"],
    sourceName = "source",
    sourceStrand =".",
    targetChrom=ccan.df[,"targetChrom"],
    targetStart=ccan.df[,"targetStart"],
    targetEnd=ccan.df[,"targetEnd"],
    targetName = "target",
    targetStrand = "+"
  )
  head(ccan.df)

# write to file in ucsc interact format for visualization in genome browser
  trackName = paste0(clusterID,"_CCAN")
  trackDescription = "kidney celltype chromatin interactions"
  comment1 = paste0("track type=interact name=", trackName, " description=\"", trackDescription, "\" color=0,0,255, interactDirectional=false maxHeightPixels=200:100:50 visibility=full")
  comment2 = paste0("browser position chr20:44,080,406-44,753,709") 

  colnames(ccan.interact)[1] <- "#chrom"
  fileConn <- file(paste0(outDir, clusterID, ".ccan.interact"), "w")
  writeLines(comment1, fileConn)
  writeLines(comment2, fileConn)
  writeLines(paste(colnames(ccan.interact), collapse="\t"), fileConn)
  write.table(ccan.interact, sep = "\t", file = fileConn, row.names = FALSE, quote = FALSE, col.names = FALSE)
  close(fileConn)
  
})  