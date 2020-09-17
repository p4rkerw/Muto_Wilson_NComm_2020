# allele specific expression in single cell data using gatk4 docker image
# the bam files are obtained from star_wasp.sh and need to be pre-filtered
# for reads that passed the WASP filter

docker pull broadinstitute/gatk:4.1.8.1

# run the container and mount the cellranger-atac bam files
docker run \
-v /g/diabneph/cellranger_rna_counts:/counts \
-v /g/reference/GRCh38-1.2.0_premrna:/ref \
-v /g/diabneph/cellranger_atac_counts/gatk/vcf:/vcf \
-v /g/diabneph/allele_spec:/outs \
-it broadinstitute/gatk:4.1.8.1

# filter out reads that did not pass WASP filter
# reads that passed the filter are tagged as vW:i:1
gatk FilterSamReads \
      -I /counts/Control_1/Control_1Aligned.sortedByCoord.out.bam \
      -O /counts/Control_1/nowasp_Control_1Aligned.sortedByCoord.out.bam \
      -T vW \
      -TV i:1 \
      --FILTER includeTagValues

# do allele specific read counting with the filtered snv and corresponding snRNA bam
gatk ASEReadCounter \
 -R ref/fasta/genome.fa \
 -I counts/Control_1/nowasp_Control_1Aligned.sortedByCoord.out.bam \
 -V vcf/no_chr.exon.snv.vcf \
 -O outs/Control_1/ase.output.table