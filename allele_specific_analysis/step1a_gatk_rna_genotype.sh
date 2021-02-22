#!/bin/bash
# this script will generate a genotype vcf using the input bam file
# to run interactively on the RIS compute1 cluster:
# export LSF_DOCKER_VOLUMES="$HOME:$HOME \
# $STORAGE1/project/cellranger_rna_counts:$HOME/counts \
# $STORAGE1/project/vcf:$HOME/vcf \
# $SCRATCH1/reference/GRCh38-2020-A.premrna/fasta:$HOME/ref \
# $STORAGE1/reference/gatk:$HOME/gatk_bundle \
# $SCRATCH1:$SCRATCH1"
# bsub -Is -G compute-parkerw -R 'rusage[mem=32GB]' -q general-interactive -a 'docker(p4rkerw/gatk:4.1.8.1)' /bin/bash

# to run detached
# export LSF_DOCKER_VOLUMES="$HOME:$HOME \
# $STORAGE1/project/cellranger_rna_counts:$HOME/counts \
# $STORAGE1/project/vcf:$HOME/vcf \
# $SCRATCH1/reference/GRCh38-2020-A.premrna/fasta:$HOME/ref \
# $STORAGE1/reference/gatk:$HOME/gatk_bundle \
# $SCRATCH1:$SCRATCH1"
# SAMPLE_ARRAY=(1 2 3 4 5)
# for SAMPLE in ${SAMPLE_ARRAY[*]}
# do
# bsub -G compute-parkerw -J "gatkc${SAMPLE}" -R 'rusage[mem=32GB]' -q general\
#  -a 'docker(broadinstitute/gatk:4.1.8.1)' -o $SCRATCH1/gatkc${SAMPLE}.out bash $REPO/allele_specific_analysis/step1_gatk_rna_genotype.sh Control_${SAMPLE}
# done

# TO RUN LOCALLY
# SCRATCH1=/g/scratch
# docker run \
# --workdir $HOME \
# -v $HOME:$HOME \
# -v g/project/cellranger_rna_counts:$HOME/counts \
# -v g/project/vcf:$HOME/vcf \
# -v g/reference/GRCh38-2020-A.premrna/fasta:$HOME/ref \
# -v g/reference/gatk:$HOME/gatk_bundle \
# -v $SCRATCH1:$SCRATCH1 \
# -e SCRATCH1="/g/scratch" \
# --rm -it p4rkerw/gatk:4.1.8.1

# make sure that miniconda is in the path
export PATH=$PATH:/gatk:/opt/miniconda/envs/gatk/bin:/gatk:/opt/miniconda/envs/gatk/bin:/opt/miniconda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

# prepare a fasta dict file using the cellranger-rna ref
if [ -a ref/genome.dict ] ; then
   echo "Sequence dictionary detected"
else
   gatk CreateSequenceDictionary -R ref/genome.fa
fi

# create an output dir
SAMPLE_NAME=$1
mkdir $SCRATCH1/gatk/
mkdir $SCRATCH1/gatk/rna
mkdir $SCRATCH1/gatk/rna/$SAMPLE_NAME
outdir=$SCRATCH1/gatk/rna/$SAMPLE_NAME

# mark duplicates
gatk MarkDuplicates \
      --I counts/version_4.0/$SAMPLE_NAME/outs/possorted_genome_bam.bam \
      --O $outdir/marked_duplicates.bam \
      --METRICS_FILE $outdir/marked_dup_metrics.txt \
      --BARCODE_TAG BC \
      --CREATE_INDEX true

# limit to specified intervals
# this will speed up analysis by eliminating alt contigs
CONTIGS=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY)      
printf '%s\n' "${CONTIGS[@]}" > /tmp/interval.list

# split reads that contain Ns in their cigar string (e.g. spanning splicing events in RNAseq data).
gatk SplitNCigarReads \
   -R ref/genome.fa \
   -I $outdir/marked_duplicates.bam \
   -L /tmp/interval.list \
   -O $outdir/cigar_marked_duplicates.bam 

# bundle files can be obtained from gatk resource bundle on google cloud
# generate base recalibration table
gatk BaseRecalibrator \
   -I $outdir/cigar_marked_duplicates.bam \
   -R ref/genome.fa \
   --known-sites gatk_bundle/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf \
   --known-sites gatk_bundle/resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.vcf.gz \
   --known-sites gatk_bundle/resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf.gz \
   --known-sites gatk_bundle/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
   -O $outdir/recal_data.table 

# apply base quality score recalibration
gatk ApplyBQSR \
   -R ref/genome.fa \
   -I $outdir/cigar_marked_duplicates.bam \
   --bqsr-recal-file $outdir/recal_data.table \
   -O $outdir/bqsr.bam 

# single sample gvcf variant calling 
# annotate the sites with dbsnp
PROCESSES=4  # Number of threads you want
for INTERVAL in ${CONTIGS[*]}; do
echo $INTERVAL
   gatk HaplotypeCaller  \
      -R ref/genome.fa \
      -I $outdir/bqsr.bam \
      -O $outdir/output.g.$INTERVAL.vcf.gz \
      -ERC GVCF \
      -D gatk_bundle/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf \
      -L $INTERVAL \
      --verbosity ERROR &
      while (( $(jobs |wc -l) >= (( ${PROCESSES} + 1 )) )); do
         sleep 0.1
      done
done

# merge vcfs from parallel haplotypecaller 
(ls $outdir/output.g.chr*.vcf.gz) > /tmp/vcf.list
gatk GatherVcfs \
   -I /tmp/vcf.list \
   -O $outdir/output.g.vcf.gz 

# index merged vcf
gatk IndexFeatureFile \
   -I $outdir/output.g.vcf.gz 

# single sample genotyping
gatk GenotypeGVCFs \
   -R ref/genome.fa \
   -V $outdir/output.g.vcf.gz \
   -O $outdir/output.vcf.gz 

# perform hard filtering using the qual-by-depth QD score and window for snp clustering
gatk VariantFiltration \
   --R ref/genome.fa \
   --V $outdir/output.vcf.gz \
   --window 35 \
   --cluster 3 \
   --filter-name "FS" \
   --filter "FS > 30.0" \
   --filter-name "QD" \
   --filter "QD < 2.0" \
   -O $outdir/$SAMPLE_NAME.filter.output.vcf.gz


