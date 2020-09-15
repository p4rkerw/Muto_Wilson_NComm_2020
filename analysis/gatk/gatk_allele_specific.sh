# allele specific expression in single cell data
# install gatk4 docker image

docker pull broadinstitute/gatk:4.1.3.0

# run the container and mount the cellranger-atac bam files
docker run \
-v /g/diabneph/cellranger_atac_counts:/counts \
-v /g/reference/refdata-cellranger-atac-GRCh38-1.2.0/fasta:/ref \
-v /g/reference/gatk:/bundle \
-it broadinstitute/gatk:4.1.3.0

# mark duplicates
gatk MarkDuplicates \
      --INPUT /counts/version_1.2/Control_1/outs/possorted_bam.bam \
      --OUTPUT /counts/gatk/marked_duplicates.bam \
      --METRICS_FILE /tmp/marked_dup_metrics.txt \
      --BARCODE_TAG BC

# prepare a fasta dict file using the cellranger-atac ref
gatk CreateSequenceDictionary -R /ref/genome.fa

# bundle files can be obtained from gatk resource bundle on google cloud
# generate base recalibration table
gatk BaseRecalibratorSpark \
   -I /counts/gatk/marked_duplicates.bam \
   -R /ref/genome.fa \
   --known-sites /bundle/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf \
   --known-sites /bundle/resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.vcf.gz \
   --known-sites /bundle/resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf.gz \
   --known-sites /bundle/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
   -O /counts/gatk/recal_data.table \
   --spark-master local[*]

# apply base quality score recalibration
gatk ApplyBQSRSpark \
   -R /ref/genome.fa \
   -I /counts/gatk/marked_duplicates.bam \
   --bqsr-recal-file /tmp/recal_data.table \
   -O /counts/gatk/output.bam \
   --spark-master local[*]

# call variants 
