# setup salmon with conda to count the fastq reads
# see https://combine-lab.github.io/salmon/getting_started/
conda create -n salmon salmon

# activate salmon 
conda activate salmon

# download the ensembl transcriptome to ~/reference
wget ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz

# build an index of the ensembl trancriptome with salmon
salmon index -t Homo_sapiens.GRCh38.cdna.all.fa.gz -i Homo_sapiens.GRCh38.cdna.all_index

# count the paired fastq files to generate a gene by count matrix
# run in directory containing paired fastq files
# use 8 threads and output to quants folder
for fn in $(ls SRR106916*_1.fastq.gz); # grab all unique fastq pairs
do
samp=$(echo $fn | sed 's/[1-2].fastq.gz//g') # grab the basename without readpair designations
echo "Processing sample ${samp}"
salmon quant -i /home/parkerw/reference/Homo_sapiens.GRCh38.cdna.all_index -l A \
         -1 ${samp}1.fastq.gz \
         -2 ${samp}2.fastq.gz \
         -p 8 --validateMappings -o quants/${samp}quant \
         --gcBias # recommended parameter in DESeq2 documentation
done 


