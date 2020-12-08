# TO RUN INTERACTIVELY:
# export LSF_DOCKER_VOLUMES="$HOME:$HOME \
# $STORAGE1/project/vcf/atac_genotype:$HOME/vcf \
# $SCRATCH1:$SCRATCH1"
# bsub -Is -G compute-parkerw -R 'rusage[mem=32GB]' -q general-interactive -a 'docker(p4rkerw/gatk:latest)' /bin/bash

# this script will process a funcotated INFO field into a csv dataframe with corresponding headers
SAMPLE_ARRAY=(1 2 3 4 5)
for SAMPLE in ${SAMPLE_ARRAY[*]}
do
	SAMPLE_NAME=Control_${SAMPLE}
	INPUT_VCF=$HOME/vcf/$SAMPLE_NAME.variants.funcotated.vcf
	OUTPUT_DIR=$HOME/vcf

	echo "Processing ${SAMPLE_NAME}"
	bcftools query -f'[%CHROM,%POS,%REF,%ALT,%GT,%FILTER\n]' $INPUT_VCF > /tmp/temp.csv
	echo "CHROM,POS,REF,ALT,GT,FILTER" | cat - /tmp/temp.csv > /tmp/anno.csv

	# this will parse a vcf annotated by funcotator and create headers
	echo "Extracting headers from INFO field"
	cat $INPUT_VCF|grep INFO|grep FUNCOTATION|cut -d: -f2|sed s'/ //g'|sed s'/">//g'|sed s'/]//g'|cut -d '|' -f 1-130|sed 's/|/,/g' > /tmp/headers.csv

	# grab the funcotation, switch to csv and and cut first 130 fields
	echo "Converting annotation in INFO field to csv"
	bcftools query -f'[%INFO/FUNCOTATION\n]' $INPUT_VCF|cut -d ',' -f 1| cut -d '|' -f 1-129|sed 's/|/,/g' > /tmp/funco.csv

	# merge to csv file
	echo "Merging headers with annotation"
	cat /tmp/headers.csv /tmp/funco.csv > /tmp/temp.csv

	echo "Generating csv annotation file"
	paste -d',' /tmp/anno.csv /tmp/temp.csv > /tmp/formatted.variant.csv

	# filter for desired annotation columns
	echo "Filtering for specified column annotations and writing file"
	awk -v cols='CHROM,POS,REF,ALT,GT,FILTER,Gencode_27_chromosome,Gencode_27_start,Gencode_27_hugoSymbol,Gencode_27_variantClassification,Gencode_27_codonChange,gnomAD_exome_AF,gnomAD_genome_AF,Gencode_27_transcriptExon' \
	'BEGIN{FS=OFS=","; nc=split(cols, a, ",")} NR==1{for (i=1; i<=NF; i++) hdr[$i]=i} \
	{for (i=1; i<=nc; i++) if (a[i] in hdr) printf "%s%s", $hdr[a[i]], (i<nc?OFS:ORS)}' /tmp/formatted.variant.csv > $OUTPUT_DIR/$SAMPLE_NAME.formatted.variant.csv

	echo "Removing temporary files"
	rm /tmp/anno.csv /tmp/temp.csv /tmp/headers.csv /tmp/funco.csv
done


