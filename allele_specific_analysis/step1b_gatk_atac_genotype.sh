#!/bin/bash
# this script will generate a genotype vcf using the input bam file
# to run interactively on the RIS compute1 cluster:
# export LSF_DOCKER_VOLUMES="$HOME:$HOME \
# $STORAGE1/project/cellranger_atac_counts:$HOME/counts \
# $STORAGE1/project/cellranger_rna_counts:$HOME/rna_counts \
# $STORAGE1/reference/refdata-cellranger-atac-GRCh38-1.2.0/fasta:$HOME/ref \
# $STORAGE1/reference/gatk:$HOME/gatk_bundle \
# $SCRATCH1:$SCRATCH1"
# bsub -Is -G compute-parkerw -R 'rusage[mem=32GB]' -q general-interactive -a 'docker(p4rkerw/gatk:4.1.8.1)' /bin/bash

# to run detached
# export LSF_DOCKER_VOLUMES="$HOME:$HOME \
# $STORAGE1/project/cellranger_atac_counts:$HOME/counts \
# $STORAGE1/project/cellranger_rna_counts:$HOME/rna_counts \
# $STORAGE1/reference/refdata-cellranger-atac-GRCh38-1.2.0/fasta:$HOME/ref \
# $STORAGE1/reference/gatk:$HOME/gatk_bundle \
# $SCRATCH1:$SCRATCH1"
# bsub -G compute-parkerw -R 'rusage[mem=32GB]' -q general -a 'docker(broadinstitute/gatk:4.1.8.1)' -o gatk2.out bash $REPO/allele_specific_analysis/step1_gatk_genotype.sh Control_2

# prepare a fasta dict file using the cellranger-atac ref
# gatk CreateSequenceDictionary -R /ref/genome.fa

# make sure that miniconda is in the path
export PATH=$PATH:/gatk:/opt/miniconda/envs/gatk/bin:/gatk:/opt/miniconda/envs/gatk/bin:/opt/miniconda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

# create an output dir
SAMPLE_NAME=$1
mkdir -p $SCRATCH1/gatk/$SAMPLE_NAME
OUTPUT_DIR=$SCRATCH1/gatk/$SAMPLE_NAME

# mark duplicates
gatk MarkDuplicates \
      --I counts/version_1.2/$SAMPLE_NAME/outs/possorted_bam.bam \
      --O $OUTPUT_DIR/marked_duplicates.bam \
      --METRICS_FILE $OUTPUT_DIR/marked_dup_metrics.txt \
      --BARCODE_TAG CB 

# bundle files can be obtained from gatk resource bundle on google cloud
# generate base recalibration table
gatk BaseRecalibrator \
   -I $OUTPUT_DIR/marked_duplicates.bam \
   -R ref/genome.fa \
   --known-sites $HOME/gatk_bundle/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf \
   --known-sites $HOME/gatk_bundle/resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.vcf.gz \
   --known-sites $HOME/gatk_bundle/resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf.gz \
   --known-sites $HOME/gatk_bundle/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
   -O $OUTPUT_DIR/recal_data.table 

# apply base quality score recalibration
gatk ApplyBQSR \
   -R ref/genome.fa \
   -I $OUTPUT_DIR/marked_duplicates.bam \
   --bqsr-recal-file $OUTPUT_DIR/recal_data.table \
   -O $OUTPUT_DIR/bqsr.bam 

# single sample gvcf variant calling 
# annotate the sites with dbsnp
gatk HaplotypeCaller  \
   -R ref/genome.fa \
   -I $OUTPUT_DIR/bqsr.bam \
   -O $OUTPUT_DIR/output.g.vcf.gz \
   -ERC GVCF \
   -D $HOME/gatk_bundle/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf

# single sample genotyping
gatk GenotypeGVCFs \
   -R ref/genome.fa \
   -V $OUTPUT_DIR/output.g.vcf.gz \
   -O $OUTPUT_DIR/output.vcf.gz

# score the variants prior to filtering
# activate the conda packages that are required for cnnscorevariants
source activate gatk
gatk CNNScoreVariants \
   -V $OUTPUT_DIR/output.vcf.gz \
   -R ref/genome.fa \
   -O $OUTPUT_DIR/annotated.vcf

# filter variants with default tranches from gatk
mkdir vcf/atac_genotype
gatk FilterVariantTranches \
   -V $OUTPUT_DIR/annotated.vcf \
   --resource gatk_bundle/resources_broad_hg38_v0_hapmap_3.3.hg38.vcf.gz \
   --resource gatk_bundle/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
   --info-key CNN_1D \
   --snp-tranche 99.95 \
   --indel-tranche 99.4 \
   -O vcf/atac_genotype/$SAMPLE_NAME.filtered.vcf  



