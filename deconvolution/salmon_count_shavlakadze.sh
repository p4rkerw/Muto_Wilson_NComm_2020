# this script will generate read counts for the bulk RNAseq from Shavlakadze et al
# see sra_explorer_fastq_download for the files obtained from SRA explorer

# # setup salmon with conda to count the fastq reads
# # see https://combine-lab.github.io/salmon/getting_started/
# conda create -n salmon salmon

# # activate salmon 
conda activate salmon

# # download the ensembl transcriptome to ~/reference
wget ftp://ftp.ensembl.org/pub/release-100/fasta/rattus_norvegicus/cdna/Rattus_norvegicus.Rnor_6.0.cdna.all.fa.gz
# # build an index of the ensembl trancriptome with salmon
salmon index -t Rattus_norvegicus.Rnor_6.0.cdna.all.fa.gz -i Rattus_norvegicus.Rnor_6.0.cdna.all_index

# count the paired fastq files to generate a gene by count matrix
# run in directory containing paired fastq files
# use 8 threads and output to quants folder
for fn in $(ls SRR8705*.fastq.gz); # grab all unique fastq pairs
do
samp=$(echo $fn | sed 's/[1-2].fastq.gz//g') # grab the basename without readpair designations
echo "Processing sample ${samp}"
salmon quant -i /mnt/g/reference/Rattus_norvegicus.Rnor_6.0.cdna.all_index -l A \
         -1 ${samp}1.fastq.gz \
         -2 ${samp}2.fastq.gz \
         -p 8 --validateMappings -o quants/${samp}quant \
         --gcBias # recommended parameter in DESeq2 documentation
done 



