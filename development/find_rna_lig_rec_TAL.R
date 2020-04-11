library(Seurat)
library(gplots)
library(mosaic)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(circlize)
library(here)

#read ligand-receptor database into R
#download if necessary at wget 
url <- "http://fantom.gsc.riken.jp/5/suppl/Ramilowski_et_al_2015/data/PairsLigRec.txt"
download.file(url, destfile = paste0("ligand_receptor/",basename(url)))
allPair <- read.table("ligand_receptor/PairsLigRec.txt", sep ="\t", header = T,quote = "")

library(Seurat)
library(harmony)
library(ggplot2)
library(dplyr)
library(Matrix)
library(BuenColors)
library(here)
set.seed(1234)

#snRNA-seq data
rnaAggr <- readRDS("cellranger_rna_prep/rnaAggr_control.rds")

#Sub-clustering
tal <- subset(rnaAggr, ident = "TAL")
tal<- FindNeighbors(tal, dims = 1:25, verbose = TRUE, reduction = "harmony")
tal <- FindClusters(tal, verbose = TRUE, resolution = 0.2, reduction = "harmony")
tal <- RunUMAP(tal, dims = 1:25, verbose = TRUE, reduction = "harmony")
DimPlot(tal, reduction = "umap", label = TRUE, pt.size = 1) + NoLegend()

#VlnPlot(tal,features = c("nCount_RNA","nFeature_RNA","percent.mt"),pt.size = 0)
#cluster1 is low-quality cells (low "nCount_RNA","nFeature_RNA" and high "percent.mt") and removed for the sub-clustering downstream analysis of TAL
tal <- subset(tal, idents = 1, invert = T)

#Re-clustering the TAL without low-quality cells
tal<- FindNeighbors(tal, dims = 1:20, verbose = TRUE, reduction = "harmony")
tal <- FindClusters(tal, verbose = TRUE, resolution = 0.2, reduction = "harmony")
tal <- RunUMAP(tal, dims = 1:20, verbose = TRUE, reduction = "harmony")

#Annotation of TAL cells
new.cluster.ids <- c("TAL1","TAL2","ATL")
names(new.cluster.ids) <- levels(tal)
tal <- RenameIdents(tal, new.cluster.ids)
tal@meta.data$subtype <- tal@active.ident

clusterID1 <- "TAL1"
clusterID2 <- "TAL2"

#subset the seurat object for glomerular cell types
sub <- subset(tal, ident = c(clusterID1, clusterID2))



#find genes expressed in designated cell types (at least 25%)
#these are not necessarily cell-type specific markers or differentially expressed
#assign the gene column to rowname because overlapping genes in different subsets have an appended .1, .2 etc.

sub_1 <- subset(sub, ident = clusterID1)
sub_2 <- subset(sub, ident = clusterID2)

marker.c1 <- rownames(sub_1)[(1-rowSums(sub_1@assays[["SCT"]]@data==0)/length(colnames(sub_1)) > 0.25)]
marker.c2 <- rownames(sub_2)[(1-rowSums(sub_2@assays[["SCT"]]@data==0)/length(colnames(sub_2)) > 0.25)]

marker.c1 <- AverageExpression(sub_1, assays = "SCT", features = marker.c1)[["SCT"]]
colnames(marker.c1) <- "avg_expression"
marker.c2 <- AverageExpression(sub_2, assays = "SCT", features = marker.c2)[["SCT"]]
colnames(marker.c2) <- "avg_expression"

#create a dataframe with all possible ligand-receptor pairs between marker1 and marker2
FindInteractions <- function(ligands, receptors) {
  #identify ligands 
  marker1Ligands <- intersect(rownames(ligands), allPair$Ligand.ApprovedSymbol)
  marker1LigMatrix <- ligands[marker1Ligands,]
  
  #identify receptors 
  marker2Receptors <- intersect(rownames(receptors), allPair$Receptor.ApprovedSymbol)
  marker2RecMatrix <- receptors[marker2Receptors,]
  
  i <- 1
  ligRecPair <- list()
  for(ligand in marker1Ligands){
    for(receptor in marker2Receptors){
      temp <- paste0(ligand,"_",receptor)
      ligRecPair[i] <- temp[1]
      i <- i + 1
    }
  }
  
  allPossiblePairs <- cbind(ligRecPair)
  
  #search the PairsLigRec.txt datafile for matches with the dataframe of possible ligand-receptor pairs
  ligRecPair <- list()
  i <- 1
  for(testPair in allPossiblePairs){
    temp <- (as.symbol(testPair)) #assign string to object
    for(pair in allPair$Pair.Name){
      if(temp == pair){
        print(testPair) #print matches
        ligRecPair[i] <- testPair
        i <- i + 1
      }
    }
  }
  ligRecPair <- unique(ligRecPair) #remove duplicates
  ligRecPair <- data.frame(cbind(ligRecPair)) #create a dataframe
  ligRecPair <- separate(data=ligRecPair,col=ligRecPair,into=c("Ligands","Receptors"),sep="_",remove=TRUE) #requires tidyr package
  
  df <- cbind(ligRecPair,ligFC=ligands[ligRecPair$Ligands,1], recFC=receptors[ligRecPair$Receptors,1])
  #df <- cbind(ligRecPair,ligLFC=5, recLFC=5)
  return(df)
}

#create matrices for all cell-type specific intercellular ligand-receptor combinations
lr_intercellular_1_to_2 <- FindInteractions(ligands=marker.c1, receptors=marker.c2)
lr_intercellular_2_to_1 <- FindInteractions(ligands=marker.c2, receptors=marker.c1)


# identify intracellular ligand-receptor interactions within individual celltypes
lr_intracellular1 <- FindInteractions(ligands=marker.c1, receptors=marker.c1)
lr_intracellular2 <- FindInteractions(ligands=marker.c2, receptors=marker.c2)


#interact <- rbind(lr_intercellular_1_to_2)
#interact<-interact[order(interact$Ligands,interact$Receptors),]


# plot the differential intracellular pathways in cluster 2
interact <- lr_intercellular_2_to_1
#exclude <- c("EGF","ERBB4")
#interact <- interact[!(interact$Ligands %in% exclude), ]
#interact <- interact[!(interact$Receptors %in% exclude), ]


circos.clear()
png(paste0("plots/ligRec.upregulated.",clusterID2,"_vs_",clusterID1,".png"))
chordDiagram(interact,
             annotationTrack = "grid",
             directional=1,
             direction.type="arrows",
             preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(interact))))))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter,
              CELL_META$ylim[1],
              CELL_META$sector.index,
              facing = "clockwise",
              niceFacing = TRUE,
              adj = c(0, 0.5))
}, bg.border = NA) # here set bg.border to NA is important
dev.off()



