#!/bin/bash
# allele specific expression in single cell data using gatk4 docker image
# the bam files are obtained from star_wasp.sh and need to be pre-filtered
# for reads that passed the WASP filter
# export LSF_DOCKER_VOLUMES="$HOME:$HOME \
# $SCRATCH1/gatk:$HOME/vcf \
# $STORAGE1/diabneph/cellranger_rna_counts:$HOME/counts \
# $STORAGE1/reference/GRCh38-1.2.0_premrna/fasta:$HOME/ref \
# $SCRATCH1:$SCRATCH1"
# bsub -Is -G compute-parkerw -R 'rusage[mem=32GB]' -q general-interactive -a 'docker(broadinstitute/gatk:4.1.8.1)' /bin/bash

#TO RUN DETACHED
# export LSF_DOCKER_VOLUMES="$HOME:$HOME \
# $SCRATCH1/gatk:$HOME/vcf \
# $STORAGE1/diabneph/cellranger_rna_counts:$HOME/counts \
# $STORAGE1/reference/GRCh38-1.2.0_premrna/fasta:$HOME/ref \
# $SCRATCH1:$SCRATCH1"
# bsub -G compute-parkerw -R 'rusage[mem=32GB]' -q general -o allspec1.out -a 'docker(broadinstitute/gatk:4.1.8.1)' bash healthy_dev/gatk/gatk_allele_specific.sh Control_1

# run the container and mount the STAR scRNA aligned bam files
# along with the vcf containing patient genotypes
# docker run \
# -v /g/diabneph/cellranger_rna_counts:/counts \
# -v /g/reference/GRCh38-1.2.0_premrna:/ref \
# -v /g/diabneph/cellranger_atac_counts/gatk/vcf:/vcf \
# -v /g/diabneph/allele_spec:/outs \
# -it broadinstitute/gatk:4.1.8.1

export PATH=$PATH:/gatk:/opt/miniconda/envs/gatk/bin:/gatk:/opt/miniconda/envs/gatk/bin:/opt/miniconda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

# prepare a fasta dict file using the cellranger ref
# gatk CreateSequenceDictionary -R /ref/fasta/genome.fa
sampleName=$1
mkdir $SCRATCH1/allele_spec
mkdir $SCRATCH1/allele_spec/$sampleName

# do allele specific read counting with the filtered snv and corresponding cellranger snRNA bam
gatk ASEReadCounter \
 -R $HOME/ref/genome.fa \
 -I $HOME/counts/version_3.1.0/$sampleName/outs/possorted_genome_bam.bam \
 -V $HOME/vcf/$sampleName/no_chr.exon.snv.vcf \
 -O $SCRATCH1/allele_spec/$sampleName/ase.output.table

# add read groups to bam files aligned with STAR to be compatible with gatk tools
# this can also be incorporated at the alignment step 
# gatk AddOrReplaceReadGroups \
# -I /counts/Control_1/nowasp.Control_1Aligned.sortedByCoord.out.bam \
# -O /counts/Control_1/rg.nowasp.Control_1Aligned.sortedByCoord.out.bam \
# -RGID Control_1 \
# -RGLB lib1 \
# -RGPL ILLUMINA \
# -RGPU Control_1 \
# -RGSM Control_1 

# perform allele specific counting
# gatk ASEReadCounter \
#  -R /ref/fasta/genome.fa \
#  -I /counts/Control_1/rg.nowasp.Control_1Aligned.sortedByCoord.out.bam \
#  -V /vcf/Control_1.vcf \
#  -O /outs/Control_1/ase.nowasp.output.table 