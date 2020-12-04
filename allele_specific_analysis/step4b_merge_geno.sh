#!/bin/bash

# TO RUN INTERACTIVE ON RIS
# export LSF_DOCKER_VOLUMES="$HOME:$HOME \
# $STORAGE1/project/vcf:$HOME/vcf \
# $SCRATCH1:$SCRATCH1"
# bsub -Is -G compute-parkerw -R 'rusage[mem=128GB]' -q general-interactive -a 'docker(p4rkerw/gatk:latest)' /bin/bash

# TO RUN DETACHED ON RIS
# export LSF_DOCKER_VOLUMES="$HOME:$HOME \
# $STORAGE1/project/vcf:$HOME/vcf \
# $SCRATCH1:$SCRATCH1"
# SAMPLE_ARRAY=(1 2 3 4 5)
# for SAMPLE in ${SAMPLE_ARRAY[*]}
# do
# bsub -G compute-parkerw -R 'rusage[mem=8GB]' -q general -a 'docker(p4rkerw/gatk:latest)' bash healthy_dev/gatk/rna_allele_specific/step4b_merge_geno.sh
# done

# TO RUN INTERACTIVE LOCALLY
# SCRATCH1=/g/scratch
# docker run \
# --workdir $HOME \
# -v $HOME:$HOME \
# -v /g/project/vcf:$HOME/vcf \
# -v $SCRATCH1:$SCRATCH1 \
# -e SCRATCH1="/g/scratch" \
# --rm -it p4rkerw/gatk:latest

# grab all the vcf snRNA and snATAC file names and put into an array
ATAC_VCF_ARRAY=$(ls vcf/atac_genotype/*vcf | awk -F'/' '{print $NF}' | grep 'mrna\|premrna' | LC_ALL=C sort)
RNA_VCF_ARRAY=$(ls vcf/rna_genotype/*vcf | awk -F'/' '{print $NF}' | grep 'mrna\|premrna' | LC_ALL=C sort)

# find intersection between the arrays
INTERSECTION=$(echo ${ATAC_VCF_ARRAY[@]} ${RNA_VCF_ARRAY[@]} | sed 's/ /\n/g' | sort | uniq -d)

# iterate over the atac vcf array and match up file names in the rna vcf array
mkdir vcf/joint_genotype
for VCF in ${INTERSECTION[*]}
do
	echo "Merging genotypes for ${VCF}"
	ATAC_VCF=vcf/atac_genotype/$VCF
	cp $ATAC_VCF /tmp/atac.vcf
	bgzip -f /tmp/atac.vcf
	tabix -f /tmp/atac.vcf.gz

	RNA_VCF=vcf/rna_genotype/$VCF
	cp $RNA_VCF /tmp/rna.vcf
	bgzip -f /tmp/rna.vcf
	tabix -f /tmp/rna.vcf.gz

	bcftools isec -p /tmp/isec /tmp/rna.vcf.gz /tmp/atac.vcf.gz > /tmp/rna_only_vcf.gz
	bgzip -f /tmp/isec/0000.vcf
	tabix -f /tmp/isec/0000.vcf.gz
	
	# merge genotypes prioritizing ATAC genotypes
	bcftools merge --force-samples /tmp/isec/0000.vcf.gz /tmp/atac.vcf.gz > vcf/joint_genotype/$VCF
	rm -rf /tmp/isec
done

