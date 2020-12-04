#!/bin/bash
# this script will annotate vcf using gatk funcotator
# to run interactively on the RIS compute1 cluster:
# export LSF_DOCKER_VOLUMES="$HOME:$HOME \
# $STORAGE1/project/cellranger_atac_counts:$HOME/counts \
# $STORAGE1/project/vcf:$HOME/vcf \
# $STORAGE1/project/cellranger_rna_counts:$HOME/rna_counts \
# $SCRATCH1/reference/GRCh38-2020-A.premrna/fasta:$HOME/ref \
# $STORAGE1/reference/gatk:$HOME/gatk_bundle \
# $SCRATCH1:$SCRATCH1"
# bsub -Is -G compute-parkerw -R 'rusage[mem=32GB]' -q general-interactive -a 'docker(p4rkerw/gatk:4.1.8.1)' /bin/bash

# to run detached
# export LSF_DOCKER_VOLUMES="$HOME:$HOME \
# $STORAGE1/project/vcf:$HOME/vcf \
# $STORAGE1/project/cellranger_atac_counts:$HOME/counts \
# $STORAGE1/project/cellranger_rna_counts:$HOME/rna_counts \
# $SCRATCH1/reference/GRCh38-2020-A.premrna/fasta:$HOME/ref \
# $STORAGE1/reference/gatk:$HOME/gatk_bundle \
# $SCRATCH1:$SCRATCH1"
# SAMPLE_ARRAY=(1 2 3)
# for SAMPLE in ${SAMPLE_ARRAY[*]}
# do
# bsub -G compute-parkerw -R 'rusage[mem=32GB]' -q general -J "annoc${SAMPLE}" -a 'docker()' -o $SCRATCH1/annoc${SAMPLE}.out \
# bash healthy_dev/gatk/rna_allele_specific/step2_gatk_anno_vcf.sh Control_${SAMPLE} $HOME/vcf/atac_genotype
# done

export PATH=$PATH:/gatk:/opt/miniconda/envs/gatk/bin:/gatk:/opt/miniconda/envs/gatk/bin:/opt/miniconda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

# create an output dir
SAMPLE_NAME=$1
VARIANT_DIR=$2

# download funcotator resource to scratch (if not already there)
# gatk FuncotatorDataSourceDownloader --germline --validate-integrity --extract-after-download

# annotate variants
gatk Funcotator \
     --variant $VARIANT_DIR/filter.output.vcf.gz \
     --reference ref/genome.fa \
     --ref-version hg38 \
     --data-sources-path $SCRATCH1/reference/funcotator_dataSources.v1.6.20190124g  \
     --output $VARIANT_DIR/$SAMPLE_NAME.variants.funcotated.vcf \
     --output-file-format VCF