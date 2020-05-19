library(Seurat) 
library(Signac) 
library(BSgenome.Hsapiens.UCSC.hg38)
library(JASPAR2018)
library(TFBSTools)
library(EnsDb.Hsapiens.v86)
library(openxlsx)
set.seed(1234)

#Use the object with all version of jasper motifs
rnaAggr <- readRDS("cellranger_rna_prep/rnaAggr_control.rds")
sub_atac <- readRDS("cellranger_atac_prep/atacAggr_sub97_control_allT.rds")
DefaultAssay(sub_atac) <- "peaks"

dac <- FindMarkers(sub_atac, ident.1 = "PT_VCAM1",
                   ident.2 = c("PCT","PST"), test.use = 'LR', latent.vars = "nCount_peaks", only.pos = TRUE)
cf <- ClosestFeature(rownames(dac), annotation=EnsDb.Hsapiens.v86, sep=c(':','-'))
dac <- cbind(dac, gene=cf$gene_name, distance=cf$distance)
dac <- dac[dac$p_val_adj < 1e-05 & dac$avg_logFC > 0.5, ]

test <- sub_atac@assays[["peaks"]]@misc[["motif"]]@data
test <- as.data.frame(test)
test$peak <- rownames(test)
test_dac <- test[rownames(dac),]
motif_interest <- c("MA0107.1","MA0105.4") #RELA/NFKB1
test_dac <- test_dac[,colnames(test_dac) %in% motif_interest]
colnames(test_dac) <- c("RELA","NFKB1")
test_dac$sum <- rowSums(test_dac)
dac <- cbind(dac,test_dac)
dac_nfkb <- dac[dac$sum > 0,]

#dac_rela.promoter$peak <- rownames(dac_rela.promoter)
#dac_nfkb$peak <- rownames(dac_nfkb)
deg <- FindMarkers(rnaAggr,ident.1 = "PT_VCAM1",ident.2 = "PT",only.pos = T)
deg$gene <- rownames(deg)
rna_diff <- deg[deg$p_val_adj<1e-5 & deg$avg_logFC > 0.5,]
atacrna_rela <- merge(dac_nfkb,rna_diff,by.x="gene",by.y="gene")
atacrna_rela <- atacrna_rela[atacrna_rela$distance<10000, ]
ab <- as.data.frame(table(as.character(atacrna_rela$gene)))
write.csv(ab,"Desktop/DNproject_results/191006_dacmotifanalysis/rela_dependent_gene_191006.csv")


#PT vs KIM1+ in control

dac_neg <- FindMarkers(sub_atac, ident.1 = c("PCT","PST"),
                   ident.2 = "PT_VCAM1", test.use = 'LR', latent.vars = "nCount_peaks", only.pos = TRUE)
cf_neg <- ClosestFeature(rownames(dac_neg), annotation=EnsDb.Hsapiens.v86, sep=c(':','-'))
dac_neg <- cbind(dac_neg, gene=cf_neg$gene_name, distance=cf_neg$distance)
dac_neg <- dac_neg[dac_neg$p_val_adj < 0.05, ]
daclist <- rownames(dac_neg)

test <- sub_atac@assays[["peaks"]]@misc[["motif"]]@data
test <- as.data.frame(test)
test$peak <- rownames(test)
test_dac <- test[daclist,]
motif_interest <- c("MA0484.1","MA0114.2","MA0046.2") #REL/RELA/NFKB1/NFKB2
test_dac <- test_dac[,colnames(test_dac) %in% motif_interest]
colnames(test_dac) <- c("HNF4G","HNF4A","HNF1A")
dac_neg <- cbind(dac_neg,test_dac)
dac_hnf4a <- dac_neg[dac_neg$HNF4A>0 & dac_neg$p_val_adj<1e-5 ,]

#dac_rela.promoter$peak <- rownames(dac_rela.promoter)
dac_hnf4a$peak <- rownames(dac_hnf4a)
deg_neg <- FindMarkers(rnaAggr,ident.1 = "PT",ident.2 = "PT_VCAM1")
deg_neg$gene <- rownames(deg_neg)
rna_diff <- deg_neg[deg_neg$avg_logFC>0.5,]
atacrna_hnf4a <- merge(dac_hnf4a,rna_diff,by.x="gene",by.y="gene")
atacrna_rela <- atacrna_rela[atacrna_rela$distance<2000, ]
ab <- as.data.frame(table(as.character(atacrna_rela$gene)))
write.csv(ab,"Desktop/DNproject_results/191006_dacmotifanalysis/rela_dependent_gene_191006.csv")













dac <- FindMarkers(atacAggr, ident.1 = "PT_control",
                   ident.2 = "PT(KIM1+)_control", test.use = 'LR', latent.vars = "nCount_peaks", only.pos = TRUE)
cf <- ClosestFeature(rownames(dac), annotation=EnsDb.Hsapiens.v86, sep=c(':','-'))
dac <- cbind(dac, gene=cf$gene_name, distance=cf$distance)
daclist <- rownames(dac)

test <- atacAggr@assays[["peaks"]]@misc[["motif"]]@data
test <- as.data.frame(test)
test$peak <- rownames(test)
test_dac <- test[daclist,]
motif_interest <- c("MA0114.1","MA0484.1","MA0046.2") #HNF4A/HNF4G/HNF1A
test_dac <- test_dac[,motif_interest]
colnames(test_dac) <- c("HNF4G","HNF1A","NR1H2::RXRA","PPARA::RXRA")
dac <- cbind(dac,test_dac)
dac_hnf4g <- dac[dac$HNF4G>0,]
#dac_rela.promoter$peak <- rownames(dac_rela.promoter)
dac_hnf4g$peak <- rownames(dac_hnf4g)
rna_diff <- read.csv("~/Desktop/DNproject_results/191005_PTvsPTKIM1/marker_cont.csv")
colnames(rna_diff) <- c("gene","p_val","avg_logFC","pct.1","pct.2","p_val_adj")
rna_diff <- rna_diff[rna_diff$avg_logFC<0,]
atacrna_hnf4g <- merge(dac_hnf4g,rna_diff,by.x="gene",by.y="gene")
atacrna_hnf4g <- atacrna_hnf4g[atacrna_hnf4g$distance<2000, ]
ab <- as.data.frame(table(as.character(atacrna_hnf4g$gene)))
write.csv(ab,"~/Desktop/DNproject_results/191006_dacmotifanalysis/hnf4g_dependent_gene_191006.csv")











