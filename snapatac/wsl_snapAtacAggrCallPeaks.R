library(Signac) # to use Extend function for GRanges objects
library(Seurat) # to perform label transfer
library(SnapATAC)
library(viridisLite)
library(ggplot2)
library(EnsDb.Hsapiens.v86) #Sys.setenv(R_INSTALL_STAGED = FALSE) to get BiocManager::install to work properly
library(GenomicRanges)
library(harmony)
library(tibble)
library(openxlsx)

x.after.sp <- readRDS(file="/mnt/g/diabneph/snapatac/snap_prep_control/atacAggr.nopeaks.snap.rds")

# identify peaks
# need to install macs2 and numpy into ubuntu which requires >= python3.6 and pip3
# cant run from RStudio in windows because it cannot access the linux executables
# call peaks for all clusters with more than 200 cells and output as a narrowpeak file
# then compile the narrowpeak files into a GRanges object peak.gr
# peak calling can be performed on clusters called by snapatac or predicted.id from snRNAseq label transfer
# the idents in the cluster slot are the predicted.id
clusters.sel = names(table(x.after.sp@cluster))[which(table(x.after.sp@cluster) > 100)]

# consider ordering the cluster list by most abundant cell type because those clusters
# will take longer to perform peak calling

peaks.ls = mclapply(seq(clusters.sel), function(i){
  print(clusters.sel[i])
  runMACS(
    obj=x.after.sp[which(x.after.sp@cluster==clusters.sel[i]),], 
    output.prefix=paste0("kidney.peaks.", gsub(" ", "_", clusters.sel)[i]),
    path.to.snaptools="/home/parkerw/anaconda3/bin/snaptools",
    path.to.macs="/home/parkerw/anaconda3/bin/macs2",
    gsize="hs", # mm, hs, etc
    buffer.size=500, 
    num.cores=4, #snaptools preprocesses each snap file separately prior to peak calling
    macs.options="--nomodel --shift 100 --ext 200 --qval 5e-2 -B --SPMR",
    tmp.folder=tempdir()
  )
}, mc.cores=3)
# assuming all .narrowPeak files in the current folder are generated from the clusters
peaks.names = system("ls | grep narrowPeak", intern=TRUE)
peak.gr.ls = lapply(peaks.names, function(x){
  peak.df = read.table(x)
  GRanges(peak.df[,1], IRanges(peak.df[,2], peak.df[,3]))
})
peak.gr = reduce(Reduce(c, peak.gr.ls))
peak.gr

# create a cell by peak matrix and add to the snap file
peaks.df = as.data.frame(peak.gr)[,1:3]
write.table(peaks.df, file = "/mnt/g/diabneph/snapatac/narrowPeak_control/peaks.combined.bed",append=FALSE,
            quote= FALSE,sep="\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"),
            fileEncoding = "")


# add the peak matrix to the snap files
system("snaptools snap-add-pmat --snap-file /mnt/g/diabneph/snapatac/snap_files/control_peaks_added/Control_1.possorted.snap --peak-file /mnt/g/diabneph/snapatac/narrowPeak_control/peaks.combined.bed")
system("snaptools snap-add-pmat --snap-file /mnt/g/diabneph/snapatac/snap_files/control_peaks_added/Control_2.possorted.snap --peak-file /mnt/g/diabneph/snapatac/narrowPeak_control/peaks.combined.bed")
system("snaptools snap-add-pmat --snap-file /mnt/g/diabneph/snapatac/snap_files/control_peaks_added/Control_3.possorted.snap --peak-file /mnt/g/diabneph/snapatac/narrowPeak_control/peaks.combined.bed")
system("snaptools snap-add-pmat --snap-file /mnt/g/diabneph/snapatac/snap_files/control_peaks_added/Control_4.possorted.snap --peak-file /mnt/g/diabneph/snapatac/narrowPeak_control/peaks.combined.bed")
system("snaptools snap-add-pmat --snap-file /mnt/g/diabneph/snapatac/snap_files/control_peaks_added/Control_5.possorted.snap --peak-file /mnt/g/diabneph/snapatac/narrowPeak_control/peaks.combined.bed")

# add cell by peak matrix to snap object
x.after.sp = readRDS("/mnt/g/diabneph/snapatac/snap_prep_control/atacAggr.nopeaks.snap.rds") # save cluster info to the snap object
x.after.sp = addPmatToSnap(x.after.sp)
x.after.sp = makeBinary(x.after.sp, mat="pmat")
x.after.sp

saveRDS(x.after.sp, file="/mnt/g/diabneph/snapatac/snap_prep_control/atacAggr.peaks.snap.rds")

# find differentially accessible regions for a single cluster compared to all
# other clusters
# the output will compare all peaks detected in all clusters and organize them
# by rows where the row number corresponds to the peak number in the peak.gr
# file calculated from runMACS. Search for genomic loci with peak.gr[peakNum]
x.after.sp@metaData$predicted.id <- as.factor(x.after.sp@metaData$predicted.id)
x.after.sp@cluster <- x.after.sp@metaData$predicted.id
dar.ls = lapply(levels(x.after.sp@cluster), function(cluster_i){
  DARs = findDAR(
    obj=x.after.sp,
    input.mat="pmat",
    cluster.pos=cluster_i,
    cluster.neg=NULL,
    cluster.neg.method="knn",
    bcv=0.4,
    test.method="exactTest",
    seed.use=10
    );
  DARs$FDR = p.adjust(DARs$PValue, method="BH");
  DARs <- rownames_to_column(DARs, var = "peak_row_index") %>%
    dplyr::filter(PValue < 0.05)
  })
names(dar.ls) = levels(x.after.sp@cluster); # idy correspond to row numbers in peak.gr

abline(h = 0, lwd=1, lty=2);
covs = Matrix::rowSums(x.after.sp@pmat);
pdf("plots/snapatac_dar.pdf")
par(mfrow = c(4,4))
for(cluster_i in levels(x.after.sp@cluster)){
  print(cluster_i)
  idy = dar.ls[[cluster_i]]$peak_row_index;
  vals = Matrix::rowSums(x.after.sp@pmat[,idy]) / covs;
  vals.zscore = (vals - mean(vals)) / sd(vals);
  plotFeatureSingle(
    obj=x.after.sp,
    feature.value=vals.zscore,
    method="umap", 
    main=cluster_i,
    point.size=0.1, 
    point.shape=19, 
    down.sample=5000,
    quantiles=c(0.01, 0.99)
    );
}
dev.off()

# modify seqnames because narrowpeak files add in leading "b"
sequpdate <- gsub("b'",'',peak.gr@seqnames) %>%
  gsub("'",'',.)
peak.gr@seqnames <- sequpdate
seqlevels(peak.gr) <- levels(peak.gr@seqnames)

# peak indices are stored in idy.ls and correspond to row numbers in peak.gr
# take a filtered list of DARs, pull genomic coords from peak.gr and annotate with gene symbol
ls.peaks.gr <- lapply(names(dar.ls), function(ident) {
  lookup <- dar.ls[[ident]]$peak_row_index
  DAR <- peak.gr[as.numeric(lookup)]
  DAR$logFC <- dar.ls[[ident]]$logFC
  DAR$PValue <- dar.ls[[ident]]$PValue
  DAR$FDR <- dar.ls[[ident]]$FDR
  DAR$logCPM <- dar.ls[[ident]]$logCPM
  DAR <- keepStandardChromosomes(DAR, pruning.mode="coarse") # filter out the alt contigs
  return(DAR)
  })
names(ls.peaks.gr) <- levels(x.after.sp@cluster)

# annotate the peaks with the nearest genomic feature
ls.anno.df <- lapply(names(ls.peaks.gr), function(ident) {
  DAR <- ls.peaks.gr[[ident]]
  chr <- as.character(DAR@seqnames)
  start <- DAR@ranges@start
  end <- DAR@ranges@start + DAR@ranges@width - 1
  lookup <- paste0(chr,":",start,"-",end)
  cf <- ClosestFeature(lookup, annotation=EnsDb.Hsapiens.v86, sep=c(':','-')) # from signac package
  DAR$gene_name <- cf$gene_name
  DAR$distance <- cf$distance
  return(as.data.frame(DAR))
  })
names(ls.anno.df) <- names(ls.peaks.gr)

# write the gr to an xlsx file
dir.create("/mnt/g/diabneph/snapatac/snap_analysis_control", showWarnings = FALSE)
write.xlsx(ls.anno.df,
 file = "/mnt/g/diabneph/snapatac/snap_analysis_control/dar.snap.celltype.control.xlsx",
 sheetName = names(ls.anno.df),
 rowNames = T)



