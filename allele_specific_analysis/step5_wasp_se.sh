#!/bin/bash
# this script will take a vcf and a cellranger bam file of single-end R2 only reads (eg. libraries 1-3)

# TO RUN INTERACTIVE ON RIS
# export LSF_DOCKER_VOLUMES="$HOME:$HOME \
# $STORAGE1/project/vcf:$HOME/vcf \
# $STORAGE1/project/cellranger_rna_counts/version_4.0:$HOME/counts \
# $STORAGE1/project/analysis/control/barcodes:$HOME/barcodes \
# $STORAGE1/reference/GRCh38-2020-A.premrna:$HOME/ref \
# $SCRATCH1:$SCRATCH1"
# bsub -Is -G compute-parkerw -R 'rusage[mem=128GB]' -q general-interactive -a 'docker(p4rkerw/wasp:1.0)' /bin/bash

# TO RUN DETACHED ON RIS
# export LSF_DOCKER_VOLUMES="$HOME:$HOME \
# $STORAGE1/project/vcf:$HOME/vcf \
# $STORAGE1/project/cellranger_rna_counts/version_4.0:$HOME/counts \
# $STORAGE1/project/analysis/control/barcodes:$HOME/barcodes \
# $STORAGE1/reference/GRCh38-2020-A.premrna:$HOME/ref \
# $SCRATCH1:$SCRATCH1"
# SAMPLE_ARRAY=(1 2 3)
# for SAMPLE in ${SAMPLE_ARRAY[*]}
# do
# bsub -G compute-parkerw -J "wasp${SAMPLE}" -o "wasp${SAMPLE}.out" -R 'rusage[mem=128GB]' -q general -a 'docker(p4rkerw/wasp:1.0)' \
# bash healthy_dev/gatk/rna_allele_specific/step5_wasp_se.sh Control_${SAMPLE} joint premrna
# done

# TO RUN INTERACTIVE LOCALLY
# SCRATCH1=/g/scratch
# docker run \
# --workdir $HOME \
# -v $HOME:$HOME \
# -v /g/project/vcf:$HOME/vcf \
# -v /g/project/analysis/control/barcodes:$HOME/barcodes \
# -v /g/project/cellranger_rna_counts/version_4.0:$HOME/counts \
# -v /g/reference/GRCh38-2020-A.premrna:$HOME/ref \
# -v $SCRATCH1:$SCRATCH1 \
# -e SCRATCH1="/g/scratch" \
# --rm -it p4rkerw/wasp:1.0 

export PATH=/opt/conda/bin:$PATH

SAMPLE_NAME=$1 # eg. Control_1
GENOTYPE=$2 # atac rna joint
VCF_REGION=$3 # mrna premrna promoter
VCF_DIR=vcf/${GENOTYPE}_genotype
CONTIG_DIR=/tmp/vcf
WASP=/WASP
OUTPUT_DIR=$SCRATCH1/wasp_snRNA_${GENOTYPE}_${VCF_REGION}
WORK_DIR=$OUTPUT_DIR/find_intersecting_snps_${SAMPLE_NAME}

mkdir $OUTPUT_DIR
mkdir $WORK_DIR
mkdir $CONTIG_DIR

# filter for sample barcodes that pass QC
# sampleName is in the 4th column
cat barcodes/barcodes.csv|grep $SAMPLE_NAME|awk -F',' '{print $1}' > /tmp/barcodes.$SAMPLE_NAME.sel.csv      

# filter the cellranger bam file by selected barcodes
BAMFILE=counts/$SAMPLE_NAME/outs/possorted_genome_bam.bam
if [ -f $OUTPUT_DIR/$SAMPLE_NAME.filter.bam ]; then
  echo "Filtered bam file detected"
else
  echo "Filtering bam file for selected barcodes"
  (samtools view -H $BAMFILE; samtools view $BAMFILE | LC_ALL=C grep -F -f /tmp/barcodes.$SAMPLE_NAME.sel.csv) | \
    samtools view -bS - > $OUTPUT_DIR/$SAMPLE_NAME.filter.bam       
fi

# generate array of unique contigs
CONTIGS=$(cat $VCF_DIR/$SAMPLE_NAME.$VCF_REGION.vcf| egrep -v "^#"| cut -f1|uniq)

# split the vcf into separate files by contig 
# print the coordinate ref and alt for each variant
for CONTIG in $CONTIGS; do
	OUTPUT_FILE=$CONTIG_DIR/$CONTIG.snps.txt.gz
	# get SNPs from VCF files:
	cat $VCF_DIR/$SAMPLE_NAME.$VCF_REGION.vcf | egrep -v "^#" |awk -v a=$CONTIG '{if($1 == a) {print $2,$4,$5}}'|gzip > $OUTPUT_FILE
done

# find reads that intersect snps
# this step may require a large amount of space depending on number of variants in vcf
echo "Finding intersecting variants and writing to ${WORK_DIR}"
python $WASP/mapping/find_intersecting_snps.py \
      --is_sorted \
      --output_dir $WORK_DIR \
      --snp_dir $CONTIG_DIR \
      $OUTPUT_DIR/$SAMPLE_NAME.filter.bam 

# if a STAR index is not detected then create one
if [ -d $SCRATCH1/starindex ]; then
  echo "STAR index detected"
else
  echo "Generating STAR index"
  STAR --runMode genomeGenerate \
  --runThreadN 20 \
  --genomeDir $SCRATCH1/starindex \
  --genomeFastaFiles ref/fasta/genome.fa \
  --sjdbGTFfile ref/genes/genes.gtf 
fi

# # use STAR to remap single end reads 
STAR --genomeDir $SCRATCH1/starindex \
 --runThreadN 20 \
 --readFilesIn <(gunzip -c $WORK_DIR/$SAMPLE_NAME.filter.remap.fq.gz) \
 --outSAMtype BAM SortedByCoordinate \
 --outFileNamePrefix $WORK_DIR/$SAMPLE_NAME.

# index the realigned bam file
samtools index -@20 $WORK_DIR/$SAMPLE_NAME.Aligned.sortedByCoord.out.bam

# filter the remapped/realigned reads and output as keep.bam
python $WASP/mapping/filter_remapped_reads.py \
      $WORK_DIR/$SAMPLE_NAME.filter.to.remap.bam \
      $WORK_DIR/$SAMPLE_NAME.Aligned.sortedByCoord.out.bam \
      $WORK_DIR/$SAMPLE_NAME.keep.bam

# merge remapped bam file
samtools merge --threads 20 $WORK_DIR/$SAMPLE_NAME.keep.merge.bam \
        $WORK_DIR/$SAMPLE_NAME.filter.keep.bam \
        $WORK_DIR/$SAMPLE_NAME.keep.bam  
samtools sort -@20 -o $OUTPUT_DIR/$SAMPLE_NAME.keep.merge.sort.bam \
         $WORK_DIR/$SAMPLE_NAME.keep.merge.bam 
samtools index $OUTPUT_DIR/$SAMPLE_NAME.keep.merge.sort.bam


# NOTE: the WASP tool does not account for UMI / barcodes when it removes duplicates
# confirm this by running samtools flagstat $OUTPUT_DIR/$SAMPLE_NAME.keep.merge.sort.bam
# and compare to samtools flagstat $OUTPUT_DIR/$SAMPLE_NAME.dedup.filter.wasp.bam
# the number of duplicates removed is greater than the number of marked duplicates
# SOLUTION: retain cellranger duplicate markings and ASEReadCounter will automatically filter them

# remove duplicates using the WASP tool which randomly selects reads for removal (instead of read with higher mapping which typically maps to ref)
# note that snRNA libraries 1-3 are R2 only whereas libraries 4-5 are paired
# python $WASP/mapping/rmdup.py $WORK_DIR/$SAMPLE_NAME.keep.merge.sort.bam $OUTPUT_DIR/$SAMPLE_NAME.dedup.filter.wasp.bam

# cleanup
# rm -rf $WORK_DIR