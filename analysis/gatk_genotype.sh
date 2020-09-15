#!/bin/bash
# this script will generate a genotype vcf using the input bam file
# to run interactively on the RIS compute1 cluster:
export LSF_DOCKER_VOLUMES="$HOME:$HOME \
$STORAGE1/diabneph/cellranger_atac_counts:$HOME/counts \
$STORAGE1/diabneph/cellranger_rna_counts:$HOME/rna_counts \
$STORAGE1/reference/refdata-cellranger-atac-GRCh38-1.2.0/fasta:$HOME/ref \
$STORAGE1/reference/gatk:$HOME/gatk_bundle \
$HOME/healthy_dev/gatk:$HOME/github_repository \
$SCRATCH1:$SCRATCH1"
bsub -Is -G compute-parkerw -R 'rusage[mem=128GB]' -q general-interactive -a 'docker(broadinstitute/gatk:4.1.8.1)' bash gatk_genotype.sh

# prepare a fasta dict file using the cellranger-atac ref
# gatk CreateSequenceDictionary -R /ref/genome.fa

# make sure that miniconda is in the path
export PATH=$PATH:/gatk:/opt/miniconda/envs/gatk/bin:/gatk:/opt/miniconda/envs/gatk/bin:/opt/miniconda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

# create an output dir
sampleName="Control_3"
mkdir $SCRATCH1/gatk/
mkdir $SCRATCH1/gatk/$sampleName
outdir=$SCRATCH1/gatk/$sampleName

# mark duplicates
gatk MarkDuplicates \
      --I counts/version_1.2/$sampleName/outs/possorted_bam.bam \
      --O $outdir/marked_duplicates.bam \
      --METRICS_FILE $outdir/marked_dup_metrics.txt \
      --BARCODE_TAG BC 

# bundle files can be obtained from gatk resource bundle on google cloud
# generate base recalibration table
gatk BaseRecalibrator \
   -I $outdir/marked_duplicates.bam \
   -R $HOME/ref/genome.fa \
   --known-sites $HOME/gatk_bundle/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf \
   --known-sites $HOME/gatk_bundle/resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.vcf.gz \
   --known-sites $HOME/gatk_bundle/resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf.gz \
   --known-sites $HOME/gatk_bundle/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
   -O $outdir/recal_data.table 

# apply base quality score recalibration
gatk ApplyBQSR \
   -R $HOME/ref/genome.fa \
   -I $outdir/marked_duplicates.bam \
   --bqsr-recal-file $outdir/recal_data.table \
   -O $outdir/bqsr.bam 

# single sample gvcf variant calling 
# annotate the sites with dbsnp
gatk HaplotypeCaller  \
   -R $HOME/ref/genome.fa \
   -I $outdir/bqsr.bam \
   -O $outdir/output.g.vcf.gz \
   -ERC GVCF \
   -D $HOME/gatk_bundle/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf

# single sample genotyping
gatk GenotypeGVCFs \
   -R $HOME/ref/genome.fa \
   -V $outdir/output.g.vcf.gz \
   -O $outdir/output.vcf.gz

# score the variants prior to filtering
# activate the conda packages that are required for cnnscorevariants
source activate gatk
gatk CNNScoreVariants \
   -V $outdir/output.vcf.gz \
   -R $HOME/ref/genome.fa \
   -O $outdir/annotated.vcf

# filter variants
gatk FilterVariantTranches \
   -V $outdir/annotated.vcf \
   --resource $HOME/gatk_bundle/resources_broad_hg38_v0_hapmap_3.3.hg38.vcf.gz \
   --resource $HOME/gatk_bundle/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
   --info-key CNN_1D \
   --snp-tranche 99 \
   --indel-tranche 99 \
   -O $outdir/filtered.vcf   

# create a file of exonic intervals with the ensembl package
# library(EnsDb.Hsapiens.v86)
# # get coord of protein coding genes 
# exon.coords <- exons(EnsDb.Hsapiens.v86, filter = ~ gene_biotype == "protein_coding", columns = "gene_name")
# seqlevelsStyle(exon.coords) <- 'UCSC'
# exon.coords <- keepStandardChromosomes(exon.coords, pruning.mode = 'coarse')
# df <- as.data.frame(exon.coords) %>%
# 	dplyr::mutate(chrom=paste0("chr",seqnames)) %>%
# 	dplyr::select(seqnames, start, end) %>%
# 	dplyr::filter(start != end)
# colnames(df)[1] <- paste0("#", colnames(df.exon[1]))
# write.table(df, row.names=FALSE, quote=FALSE, file="G:/downloads/hg38.noalt.exon.intervals.bed")

gatk SelectVariants \
	-R $HOME/ref/genome.fa \
	-V $outdir/filtered.vcf \
	--select-type-to-include SNP \
	-L $HOME/gatk_bundle/hg38.noalt.exon.intervals.bed \
	-O $outdir/exon.snv.vcf   

# prepare a fasta dict file using the cellranger ref
# gatk CreateSequenceDictionary -R $HOME/rna_ref/fasta/genome.fa

# change the contig style from hg38 to grch38 by removing "chr"
awk '{gsub(/^chr/,""); print}' $outdir/exon.snv.vcf | sed '/^##contig/d' > $outdir/no_chr.exon.snv.vcf

# index the updated vcf file
gatk IndexFeatureFile -I $outdir/no_chr.exon.snv.vcf
