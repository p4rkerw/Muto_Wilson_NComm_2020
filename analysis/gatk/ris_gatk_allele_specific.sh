# to run interactively on the RIS compute1 cluster:
export LSF_DOCKER_VOLUMES="$HOME:$HOME \
$STORAGE1/diabneph/cellranger_atac_counts:$HOME/counts \
$STORAGE1/diabneph/cellranger_rna_counts:$HOME/rna_counts \
$STORAGE1/reference/refdata-cellranger-atac-GRCh38-1.2.0/fasta:$HOME/ref \
$STORAGE1/reference/GRCh38-1.2.0_premrna:$HOME/rna_ref \
$STORAGE1/reference/gatk:$HOME/gatk_bundle \
$SCRATCH1:$SCRATCH1"
bsub -Is -G compute-parkerw -R 'rusage[mem=128GB]' -q general-interactive -a 'docker(broadinstitute/gatk:4.1.8.1)' /bin/bash 

# prepare a fasta dict file using the cellranger-atac ref
# gatk CreateSequenceDictionary -R /ref/genome.fa

# make sure that miniconda is in the path
export PATH=$PATH:/gatk:/opt/miniconda/envs/gatk/bin:/gatk:/opt/miniconda/envs/gatk/bin:/opt/miniconda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

# mark duplicates
gatk MarkDuplicates \
      --I counts/version_1.2/Control_1/outs/possorted_bam.bam \
      --O $SCRATCH1/marked_duplicates.bam \
      --METRICS_FILE $SCRATCH1/marked_dup_metrics.txt \
      --BARCODE_TAG BC 

# bundle files can be obtained from gatk resource bundle on google cloud
# generate base recalibration table
gatk BaseRecalibrator \
   -I $SCRATCH1/marked_duplicates.bam \
   -R $HOME/ref/genome.fa \
   --known-sites $HOME/gatk_bundle/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf \
   --known-sites $HOME/gatk_bundle/resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.vcf.gz \
   --known-sites $HOME/gatk_bundle/resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf.gz \
   --known-sites $HOME/gatk_bundle/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
   -O $SCRATCH1/recal_data.table 

# apply base quality score recalibration
gatk ApplyBQSR \
   -R $HOME/ref/genome.fa \
   -I $SCRATCH1/marked_duplicates.bam \
   --bqsr-recal-file $SCRATCH1/recal_data.table \
   -O $SCRATCH1/bqsr.bam 

# single sample gvcf variant calling 
# annotate the sites with dbsnp
gatk HaplotypeCaller  \
   -R $HOME/ref/genome.fa \
   -I $SCRATCH1/bqsr.bam \
   -O $SCRATCH1/output.g.vcf.gz \
   -ERC GVCF \
   -D $HOME/gatk_bundle/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf

# single sample genotyping
gatk GenotypeGVCFs \
   -R $HOME/ref/genome.fa \
   -V $SCRATCH1/output.g.vcf.gz \
   -O $SCRATCH1/output.vcf.gz

# score the variants prior to filtering
# activate the conda packages that are required for cnnscorevariants
source activate gatk
gatk CNNScoreVariants \
   -V $SCRATCH1/output.vcf.gz \
   -R $HOME/ref/genome.fa \
   -O $SCRATCH1/annotated.vcf

# filter variants
gatk FilterVariantTranches \
   -V $SCRATCH1/annotated.vcf \
   --resource $HOME/gatk_bundle/resources_broad_hg38_v0_hapmap_3.3.hg38.vcf.gz \
   --resource $HOME/gatk_bundle/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
   --info-key CNN_1D \
   --snp-tranche 99 \
   --indel-tranche 99 \
   -O $SCRATCH1/filtered.vcf   

# select exonic snvs 
# use intervals downloaded from ensembl gtf for hg38
# filter out the alt contigs
# df <- fread("downloads/Homo_sapiens.GRCh38.100.gtf.gz")
# colnames(df)[1] <- "chrom"
# contigs <- c(seq(1:22), "X","Y")
# df.exon <- dplyr::filter(df, grepl("protein_coding", V9)) %>%
# dplyr::filter(V3 == "exon") %>%
# dplyr::filter(chrom %in% contigs) %>%
# dplyr::filter(V2 == "ensembl") %>%
# dplyr::mutate(chrom=paste0("chr",chrom), start=V4, end=V5) %>%
# dplyr::filter(end > start) %>% 
# dplyr::select(chrom, start, end)
# colnames(df.exon)[1] <- paste0("#", colnames(df.exon[1]))
# write.table(df.exon, row.names=FALSE, quote=FALSE, file="downloads/hg38.noalt.exon.intervals.bed")
gatk SelectVariants \
	-R $HOME/ref/genome.fa \
	-V $SCRATCH1/filtered.vcf \
	--select-type-to-include SNP \
	-L $HOME/gatk_bundle/hg38.noalt.exon.intervals.bed \
	-O $SCRATCH1/exon.snv.vcf   

# prepare a fasta dict file using the cellranger ref
# gatk CreateSequenceDictionary -R $HOME/rna_ref/fasta/genome.fa

# change the contig style from hg38 to grch38 by removing "chr"
awk '{gsub(/^chr/,""); print}' $SCRATCH1/exon.snv.vcf | sed '/^##contig/d' > $SCRATCH1/no_chr.exon.snv.vcf

# index the updated vcf file
gatk IndexFeatureFile -I $SCRATCH1/no_chr.exon.snv.vcf

# do allele specific read counting with the filtered snv and corresponding snRNA bam
gatk ASEReadCounter \
 -R $HOME/rna_ref/fasta/genome.fa \
 -I $HOME/rna_counts/version_3.1.0/Cont1/outs/possorted_genome_bam.bam \
 -V $SCRATCH1/no_chr.exon.snv.vcf \
 -O $SCRATCH1/ase.output.table

