#!/bin/bash
# allele specific expression in single cell data using gatk4 docker image
# bam files are obtained from the wasp pipeline

# TO RUN INTERACTIVE ON RIS
# export LSF_DOCKER_VOLUMES="$HOME:$HOME \
# $STORAGE1/project/vcf:$HOME/vcf \
# $STORAGE1/project/analysis/control/barcodes:$HOME/barcodes \
# $STORAGE1/project/analysis/control:$HOME/outs \
# $STORAGE1/reference/GRCh38-2020-A.premrna:$HOME/ref \
# $STORAGE1/project/cellranger_rna_counts/version_4.0:$HOME/counts \
# $SCRATCH1:$SCRATCH1"
# bsub -Is -G compute-parkerw -R 'rusage[mem=32GB]' -q general-interactive -a 'docker(broadinstitute/gatk:4.1.8.1)' /bin/bash

#TO RUN DETACHED ON RIS
# export LSF_DOCKER_VOLUMES="$HOME:$HOME \
# $STORAGE1/project/vcf:$HOME/vcf \
# $STORAGE1/project/analysis/control/barcodes:$HOME/barcodes \
# $STORAGE1/reference/GRCh38-2020-A.premrna:$HOME/ref \
# $STORAGE1/project/cellranger_rna_counts/version_4.0:$HOME/counts \
# $SCRATCH1:$SCRATCH1"
# SAMPLE_ARRAY=(1 2 3 4 5)
# for SAMPLE in ${SAMPLE_ARRAY[*]}
# do
# bsub -G compute-parkerw -R 'rusage[mem=64GB]' -J allCnt_${SAMPLE} -q general -o allcnt_${SAMPLE}.out -a 'docker(broadinstitute/gatk:4.1.8.1)' \
# bash healthy_dev/gatk/rna_allele_specific/step6_gatk_alleleCount.sh Control_${SAMPLE} joint
# done

# TO RUN INTERACTIVE LOCALLY
# SCRATCH1=/g/scratch
# docker run \
# --workdir $HOME \
# -v $HOME:$HOME \
# -v /g/reference/GRCh38-2020-A.premrna:$HOME/ref \
# -v /g/project/vcf:$HOME/vcf \
# -v /g/project/analysis/control:$HOME/outs \
# -v /g/project/analysis/control/barcodes:$HOME/barcodes \
# -v $SCRATCH1:$SCRATCH1 \
# -e SCRATCH1="/g/scratch" \
# --rm -it broadinstitute/gatk:4.1.8.1

# update path for gatk and gatk conda dependencies (eg. CNN)
export PATH=$PATH:/gatk:/opt/miniconda/envs/gatk/bin:/gatk:/opt/miniconda/envs/gatk/bin:/opt/miniconda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

# prepare a fasta dict file using the cellranger ref
# gatk CreateSequenceDictionary -R /ref/fasta/genome.fa
SAMPLE_NAME=$1
GENOTYPE=$2
VCF_DIR=vcf/${GENOTYPE}_genotype
INPUT_DIR=$SCRATCH1/wasp_snRNA_joint_premrna
OUTPUT_DIR=$SCRATCH1/wasp_snRNA_joint_premrna/allele_rna_count
BAMFILE=$INPUT_DIR/$SAMPLE_NAME.keep.merge.sort.bam
mkdir -p $OUTPUT_DIR/$SAMPLE_NAME

# select SNPs (this is the only variant type accepted by ASEReadCounter)
# exclude multiallelic calls with restrict-alleles-to BIALLELIC
gatk SelectVariants \
	-R ref/fasta/genome.fa \
	-V $VCF_DIR/$SAMPLE_NAME.mrna.vcf \
	--select-type-to-include SNP \
	--restrict-alleles-to BIALLELIC \
	-O /tmp/$SAMPLE_NAME.snp.mrna.vcf
bgzip /tmp/$SAMPLE_NAME.snp.mrna.vcf
tabix /tmp/$SAMPLE_NAME.snp.mrna.vcf.gz	

gatk SelectVariants \
	-R ref/fasta/genome.fa \
	-V $VCF_DIR/$SAMPLE_NAME.premrna.vcf \
	--select-type-to-include SNP \
	--restrict-alleles-to BIALLELIC \
	-O /tmp/$SAMPLE_NAME.snp.premrna.vcf
bgzip /tmp/$SAMPLE_NAME.snp.premrna.vcf
tabix /tmp/$SAMPLE_NAME.snp.premrna.vcf.gz		

# perform allele specific counting for all celltypes (pseudobulk)
# tool will automatically filter out non-heterozygous positions and duplicates

gatk ASEReadCounter \
 -R ref/fasta/genome.fa \
 -I $BAMFILE  \
 -V /tmp/$SAMPLE_NAME.snp.mrna.vcf.gz \
 -O $OUTPUT_DIR/$SAMPLE_NAME/ase.mrna.ALLCELLS.output.table 

gatk ASEReadCounter \
 -R ref/fasta/genome.fa \
 -I $BAMFILE  \
 -V /tmp/$SAMPLE_NAME.snp.premrna.vcf.gz \
 -O $OUTPUT_DIR/$SAMPLE_NAME/ase.premrna.ALLCELLS.output.table 

# create another bam file that only contains reads for a specified cell type
# read in the barcodes with the following format:barcode,celltype,lowres.celltype,orig.ident
# ase counting for individual celltypes
CELLTYPE_GROUPS=$(tail -n+2 barcodes/barcodes.csv|awk -F',' '{print $3}'|sort|uniq)
for CELLTYPE in $CELLTYPE_GROUPS
do
	# filter barcodes.csv for selected celltype barcodes
	echo $CELLTYPE
	awk -v a="${SAMPLE_NAME}" -F',' '$4 == a' barcodes/barcodes.csv | \
	awk -v a="${CELLTYPE}" -F',' '$3 == a' | \
	awk -F',' '{print $1}' > /tmp/barcodes.$CELLTYPE.sel.csv
	
	# filter bam file for selected celltype barcodes
	echo "Filtering bam file for celltype barcodes" 
	(samtools view -H $BAMFILE; samtools view $BAMFILE | LC_ALL=C grep -F -f /tmp/barcodes.$CELLTYPE.sel.csv) | \
	samtools view -bS - > $INPUT_DIR/$SAMPLE_NAME.$CELLTYPE.keep.merge.sort.bam

	# perform celltype allele specific counting
	# tool will automatically filter out non-heterozygous positions
	gatk ASEReadCounter \
	 -R ref/fasta/genome.fa \
	 -I $INPUT_DIR/$SAMPLE_NAME.$CELLTYPE.keep.merge.sort.bam \
	 -V /tmp/$SAMPLE_NAME.snp.mrna.vcf.gz \
	 -O $OUTPUT_DIR/$SAMPLE_NAME/ase.mrna.$CELLTYPE.output.table   

	gatk ASEReadCounter \
	 -R ref/fasta/genome.fa \
	 -I $INPUT_DIR/$SAMPLE_NAME.$CELLTYPE.keep.merge.sort.bam \
	 -V /tmp/$SAMPLE_NAME.snp.premrna.vcf.gz \
	 -O $OUTPUT_DIR/$SAMPLE_NAME/ase.premrna.$CELLTYPE.output.table
done

# copy files to storage
mkdir $HOME/outs/wasp_rna
cp -r $OUTPUT_DIR/ $HOME/outs/wasp_rna
