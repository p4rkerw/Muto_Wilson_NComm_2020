library(Signac) # to use Extend function for GRanges objects
library(Seurat) # to perform label transfer
library(SnapATAC)
library(viridisLite)
library(ggplot2)
library(EnsDb.Hsapiens.v86) #Sys.setenv(R_INSTALL_STAGED = FALSE) to get BiocManager::install to work properly
library(GenomicRanges)
library(harmony)
library(tibble)
library(dplyr)
library(stringr)

# create the snap file from preprocessed 10X snATACseq output (see snapAtacPreprocess.R)
file.list = c("Control_1.possorted.snap",
              "Control_2.possorted.snap",
              "Control_3.possorted.snap",
              "Control_4.possorted.snap",
              "Control_5.possorted.snap"
              )

file.path = paste0("/mnt/g/diabneph/snapatac/snap_files/", file.list)

sample.list = c("Control_1",
                "Control_2",
                "Control_3",
                "Control_4",
                "Control_5"
                )

# create a list of snap objects to aggregate
x.sp.ls = lapply(seq(file.path), function(i){
  x.sp = createSnap(file=file.path[i], sample=sample.list[i]);
  x.sp
})
names(x.sp.ls) = sample.list;
sample.list

# extract barcodes from filtered aggregate snATAC-seq object
# barcodes have the format: barcode-1, barcode-2 etc with unique sample ID corresponding to the number
# they are stored in the following location in the seurat object:
# atacAggr@assays$peaks@counts@Dimnames[[2]]
# and have been saved as barcodes_retained_control.csv
barcodes_retained <- read.csv("snapatac/snap_prep_control/barcodes_retained_control.csv")

# filter the barcodes using the retained list
x.sp.list = lapply(seq(x.sp.ls), function(i){
  x.sp = x.sp.ls[[i]]
  barcodes <- dplyr::filter(barcodes_retained, grepl(i,x)) # grab the barcodes with the suffix corresponding to index eg. -1, -2 etc
  barcodes <- str_replace_all(barcodes$x, pattern="[0-9]", "1") # switch index to 1 so it matches with barcodes from invidividual libs
  x.sp = x.sp[which(x.sp@barcode %in% barcodes),] # filter libs for corresponding retained barcodes
})
names(x.sp.list) = sample.list

# add cell by bin matrix
x.sp.list = lapply(seq(x.sp.list), function(i){
  x.sp = addBmatToSnap(x.sp.list[[i]], bin.size=5000); # cannot select bin.size=1000 due to memory constraint
  x.sp
})
x.sp.list

# combine snap objects
bin.shared = Reduce(intersect, lapply(x.sp.list, function(x.sp) x.sp@feature$name))
x.sp.list <- lapply(x.sp.list, function(x.sp){
  idy = match(bin.shared, x.sp@feature$name)
  x.sp[,idy, mat="bmat"]
})
x.sp = Reduce(snapRbind, x.sp.list)
rm(x.sp.list); # free memory
gc()
table(x.sp@sample)

# binarize matrix
x.sp = makeBinary(x.sp, mat="bmat")

# filter bins
# system("wget http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/hg38-human/hg38.blacklist.bed.gz")
black_list = read.table("snapatac/hg38.blacklist.bed.gz")
black_list.gr = GRanges(
  black_list[,1], 
  IRanges(black_list[,2], black_list[,3])
)
idy = queryHits(findOverlaps(x.sp@feature, black_list.gr))
if(length(idy) > 0){x.sp = x.sp[,-idy, mat="bmat"]}
x.sp

# remove unwanted chr
chr.exclude = seqlevels(x.sp@feature)[grep("random|chrM", seqlevels(x.sp@feature))]
idy = grep(paste(chr.exclude, collapse="|"), x.sp@feature)
if(length(idy) > 0){x.sp = x.sp[,-idy, mat="bmat"]}
x.sp

# remove top 5% of bins
bin.cov = log10(Matrix::colSums(x.sp@bmat)+1);
bin.cutoff = quantile(bin.cov[bin.cov > 0], 0.95);
idy = which(bin.cov <= bin.cutoff & bin.cov > 0);
x.sp = x.sp[, idy, mat="bmat"];
x.sp

# reduce dimensionality
row.covs = log10(Matrix::rowSums(x.sp@bmat)+1)
row.covs.dens = density(
  x = row.covs, 
  bw = 'nrd', adjust = 1
)
sampling_prob = 1 / (approx(x = row.covs.dens$x, y = row.covs.dens$y, xout = row.covs)$y + .Machine$double.eps) 
set.seed(1)
idx.landmark.ds = sort(sample(x = seq(nrow(x.sp)), size = 10000, prob = sampling_prob))
x.landmark.sp = x.sp[idx.landmark.ds,]
x.query.sp = x.sp[-idx.landmark.ds,]
x.landmark.sp = runDiffusionMaps(
  obj= x.landmark.sp,
  input.mat="bmat", 
  num.eigs=50
)
x.query.sp = runDiffusionMapsExtension(
  obj1=x.landmark.sp, 
  obj2=x.query.sp,
  input.mat="bmat"
)
x.landmark.sp@metaData$landmark = 1
x.query.sp@metaData$landmark = 0
x.sp = snapRbind(x.landmark.sp, x.query.sp)
## combine landmarks and query cells;
x.sp = x.sp[order(x.sp@sample),] # IMPORTANT
rm(x.landmark.sp, x.query.sp) # free memory

# determine number of dimensions
pdf("plots/aggr.dimReductPW.pdf")
plotDimReductPW(
  obj=x.sp, 
  eigs.dims=1:50,
  point.size=0.3,
  point.color="grey",
  point.shape=19,
  point.alpha=0.6,
  down.sample=5000,
  pdf.file.name=NULL, 
  pdf.height=7, 
  pdf.width=7
)
dev.off()

# remove batch effect with harmony
x.after.sp = runHarmony(
  obj=x.sp, 
  eigs.dim=1:21, 
  meta_data=x.sp@sample # sample index
)

# perform graph based clustering and visualization
x.after.sp = runKNN(
  obj= x.after.sp,
  eigs.dim=1:21,
  k=20 # seurat uses k=20
)
x.after.sp = runCluster(
  obj=x.after.sp,
  tmp.folder=tempdir(),
  louvain.lib="R-igraph",
  path.to.snaptools=NULL,
  seed.use=10,
  resolution = 0.1
)
x.after.sp@metaData$cluster = x.after.sp@cluster


# pass cosine distance matrix to umap to avoid outliers. This is the standard method in the RunUMAP function
# in the seurat package
# x.sp = runViz(
#   obj=x.sp, 
#   tmp.folder=tempdir(),
#   dims=2,
#   eigs.dims=1:20, 
#   method="umap",
#   seed.use=10,
#   num.cores=10,
#   metric="cosine"
# )
x.after.sp = runViz(
  obj=x.after.sp, 
  tmp.folder=tempdir(),
  dims=2,
  eigs.dims=1:21, 
  method="umap", 
  seed.use=10,
  num.cores=10,
  metric="cosine"
)

# draw plots for before and after harmony batch effect correction
pdf("plots/harmony.aggr.umap.pdf")
par(mfrow = c(2, 3));
# plotViz(
#   obj=x.sp,
#   method="umap", 
#   main="Before Harmony",
#   point.color=x.sp@sample, 
#   point.size=0.1, 
#   text.add= FALSE,
#   down.sample=10000,
#   legend.add=TRUE
# )
plotViz(
  obj=x.after.sp,
  method="umap", 
  main="After Harmony",
  point.color=x.sp@sample, 
  point.size=0.1, 
  text.add=FALSE,
  down.sample=10000,
  legend.add=TRUE
)
plotViz(
  obj=x.after.sp,
  method="umap", 
  main="Cluster",
  point.color=x.after.sp@cluster, 
  point.size=0.1, 
  text.add=FALSE,
  text.size=1,
  text.color="black",
  text.halo.add=TRUE,
  text.halo.color="white",
  text.halo.width=0.2,
  down.sample=10000,
  legend.add=FALSE
)
dev.off()


# create GRanges object for protein coding and promoter coords
gene.coords <- genes(EnsDb.Hsapiens.v86, filter = ~ gene_biotype == "protein_coding")
seqlevelsStyle(gene.coords) <- 'UCSC'
genebody.coords <- keepStandardChromosomes(gene.coords, pruning.mode = 'coarse')
genes.gr <- Extend(x = gene.coords, upstream = 2000, downstream = 0)

# select marker genes for kidney
marker.genes <- c("CUBN","HAVCR1","SLC5A1","SLC5A2", # PT and PT-KIM1+ markers
                      "CFH", # PEC
                      "SLC12A1", # TAL NKCC2
                      "SLC12A3","TRPM6", # DCT1 and DCT2 NCC
                      "SCNN1G","TRPV5", # DCT2/CNT ENaC
                      "CALB1", # CNT
                      "AQP2", # PC
                      "ATP6V0D2", # ICA and ICB
                      "SLC4A1","SLC26A7", # ICA
                      "SLC26A4", # ICB
                      "NPHS1","NPHS2", # PODO
                      "PECAM1","FLT1","EMCN", # ENDO
                      "CLDN5", # GEC
                      "ITGA8","PDGFRB", # MES
                      "ACTA2","CALD1", # FIB
                      "PTPRC") # WBC

# overlap the marker genes
genes.sel.gr <- genes.gr[which(genes.gr$gene_name %in% marker.genes)]
names(mcols(genes.sel.gr))[2] <- "name"

# re-add the cell-by-bin matrix to the snap object
x.after.sp = addBmatToSnap(x.after.sp, do.par=TRUE, num.cores=2)
x.after.sp = createGmatFromMat(
  obj=x.after.sp, 
  input.mat="bmat",
  genes=genes.sel.gr,
  do.par=TRUE,
  num.cores=10
)

# normalize the cell-by-gene matrix
x.after.sp = scaleCountMatrix(
  obj=x.after.sp, 
  cov=log10(Matrix::rowSums(x.after.sp@bmat)+1),
  mat="gmat",
  method = "RPM"
)

# smooth the cell-by-gene matrix
# this function can cause a segfault for large datasets. See alternative function myRunMagic below
# x.after.sp = runMagic(
#   obj=x.after.sp,
#   input.mat="gmat",
#   step.size=3
# )

myRunMagic <- function (obj, input.mat, step.size) {
    A = obj@graph@mat;
    data.use = obj@gmat;
    
    # smooth
    A = A + t(A);
    A = A / Matrix::rowSums(A);
    data.use.smooth = A %*% data.use;
    if(step.size > 1){
        for(i in 1:step.size){
            data.use.smooth = A %*% data.use.smooth;
        }
    }
    
    slot(obj, input.mat) = data.use.smooth;    
    return(obj)
}

x.after.sp = myRunMagic(x.after.sp, input.mat = "gmat", step.size = 3)

pdf("plots/harmony.markerGenes.pdf")
par(mfrow = c(3, 2))
for(i in 1:length(marker.genes)){
  plotFeatureSingle(
    obj=x.after.sp,
    feature.value=x.after.sp@gmat[, marker.genes[i]],
    method="umap", 
    main=marker.genes[i],
    point.size=0.1, 
    point.shape=19, 
    down.sample=10000,
    quantiles=c(0, 1)
  )}
dev.off()

# calculate the ensemble signals for each cluster
ensemble.ls = lapply(split(seq(length(x.after.sp@cluster)), x.after.sp@cluster), function(x){
  SnapATAC::colMeans(x.after.sp[x,], mat="bmat")
})
# cluster using 1-cor as distance  
hc = hclust(as.dist(1 - cor(t(do.call(rbind, ensemble.ls)))), method="ward.D2")
pdf("plots/harmony.umap_kidney.pdf")
plotViz(
  obj=x.after.sp,
  method="umap", 
  main="10X Kidney",
  point.color=x.after.sp@cluster, 
  point.size=0.1, 
  point.shape=19, 
  point.alpha=0.8, 
  text.add=TRUE,
  text.size=1.5,
  text.color="black",
  text.halo.add=TRUE,
  text.halo.color="white",
  text.halo.width=0.2,
  down.sample=10000,
  legend.add=FALSE
)
plot(hc, hang=-1, xlab="")
dev.off()

# perform label transfer from previously annotated snRNAseq object
rnaAggr = readRDS("/mnt/g/diabneph/cellranger_rna_prep/rnaAggr_control.rds")
rnaAggr$tech = "rna"
variable.genes = VariableFeatures(object = rnaAggr)

# obtain gencode annotation from ucsc genome 
names(mcols(genes.gr))[2] <- "name"
genes.sel.gr = genes.gr[which(genes.gr$name %in% variable.genes)];
## reload the bmat, this is optional but highly recommended
x.after.sp = addBmatToSnap(x.after.sp) #, do.par=TRUE, num.cores=10);
x.after.sp = createGmatFromMat(
  obj=x.after.sp, 
  input.mat="bmat",
  genes=genes.sel.gr,
  do.par=TRUE,
  num.cores=10
)

# convert snap to seurat object and integrate
x.after.sp_seurat <- snapToSeurat(
  obj=x.after.sp, 
  eigs.dims=1:21, 
  norm=TRUE,
  scale=TRUE
)
transfer.anchors <- FindTransferAnchors(
  reference = rnaAggr, 
  query = x.after.sp_seurat, 
  features = variable.genes, 
  reference.assay = "RNA", 
  query.assay = "ACTIVITY", 
  reduction = "cca"
)
celltype.predictions <- TransferData(
  anchorset = transfer.anchors, 
  refdata = rnaAggr@active.ident,
  weight.reduction = x.after.sp_seurat[["SnapATAC"]],
  dims = 1:21
)
x.after.sp@metaData$predicted.id = celltype.predictions$predicted.id
x.after.sp@metaData$predict.max.score = apply(celltype.predictions[,-1], 1, max)
x.after.sp@cluster = as.factor(x.after.sp@metaData$predicted.id)

pdf("plots/harmony.predicted_id.pdf")
plotViz(
  obj=x.after.sp,
  method="umap", 
  main="10X Kidney",
  point.color=x.after.sp@cluster, 
  point.size=0.1, 
  point.shape=19, 
  point.alpha=0.8, 
  text.add=TRUE,
  text.size=1.5,
  text.color="black",
  text.halo.add=TRUE,
  text.halo.color="white",
  text.halo.width=0.2,
  down.sample=10000,
  legend.add=FALSE
)
dev.off()

saveRDS(x.after.sp, file="/mnt/g/diabneph/snapatac/snap_prep_control/atacAggr.nopeaks.snap.rds")
