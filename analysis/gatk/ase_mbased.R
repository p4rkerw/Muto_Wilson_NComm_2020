library(MBASED)
library(data.table)
library(dplyr)
library(stringr)
library(GenomicRanges)
library(EnsDb.Hsapiens.v86)
library(VariantAnnotation)

# read in the allele-specific expression read counts from snRNA ASEReadCounter following GATK identification
# of biallelic SNVs from snATAC data
ase <- fread("cellranger_atac_counts/gatk/ase.output.table")

# load the snv from the count table into a GRanges object
snv <- GRanges(
  seqnames=paste0("chr",ase$contig),
  ranges=IRanges(start=ase$position, width=1),
  allele1=ase$refAllele,
  allele2=ase$altAllele,
  allele1_counts=ase$refCount,
  allele2_counts=ase$altCount
)

# load the ensembl annotations
gene.coords <- genes(EnsDb.Hsapiens.v86, filter = ~ gene_biotype == "protein_coding", columns = "gene_name")
seqlevelsStyle(gene.coords) <- 'UCSC'
gene.coords <- keepStandardChromosomes(gene.coords, pruning.mode = 'coarse')

# find the genes that correspond to each snv (some snv overlap with multiple genes)
snv <- join_overlap_intersect(snv, gene.coords)
snv$aseID <- snv$gene_name

# create a summarized experiment 
se <- SummarizedExperiment(
  assays=list(
    lociAllele1Counts=matrix(
      snv$allele1_counts,
      ncol=1,
      dimnames=list(
        seq(1:length(snv)), 
        'project'
      )
    ),
    lociAllele2Counts=matrix(
      snv$allele2_counts,
      ncol=1, 
      dimnames=list(
        seq(1:length(snv)), 
        'project'
      )
    )
  ),
  rowRanges=snv
)


# use MBASED to evaluate allele specific expression
ASEresults_1s_haplotypesKnown <- runMBASED(
  ASESummarizedExperiment=se,
  isPhased=FALSE,
  numSim=10^6,
  BPPARAM = SerialParam()
)

res = assays(ASEresults_1s_haplotypesKnown)
majorAF <- as.data.frame(res$majorAlleleFrequency)
pval <- data.frame(res$pValueASE)
phet <- pval <- data.frame(res$pValueHeterogeneity)
df <- data.frame(gene=rownames(majorAF), majorAF=majorAF$project, pval=pval$project, phet=phet$project)

# adjust the pvalues
library(multtest)
padj <- mt.rawp2adjp(df$pval, proc = "BH")
padj <- padj$adjp[order(padj$index),]
df$padj <- padj[,2] # this is the benjamini hochberg adjusted pval

df.sig <- dplyr::filter(df, padj < 0.05)


