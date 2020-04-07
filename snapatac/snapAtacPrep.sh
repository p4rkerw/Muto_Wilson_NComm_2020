# this script will convert 10x genomics snATACseq output into a snap file that is
# compatible with the SnapATAC workflow
# to run:
# snapAtacPrep.sh LibraryName
# eg. bash snapAtacPrep.sh Control_1
# extract header file
samtools view /mnt/g/scratch/$1/outs/possorted_bam.bam -H > $1_possorted.header.sam

cat <( cat $1_possorted.header.sam) \
<( samtools view /mnt/g/scratch/$1/outs/possorted_bam.bam \
|awk '{for (i=12; i<=NF; ++i) { if ($i ~ "^CB:Z:"){ td[substr($i,1,2)] = substr($i,6,length($i)-5); } }; printf "%s:%s\n", td["CB"], $0 }' ) \
| samtools view -bS - > $1.possorted.snap.bam

samtools sort -n -@ 10 -m 4G $1.possorted.snap.bam -o $1.possorted.snap.nsrt.bam

# wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes

# generate snap file
snaptools snap-pre  \
	--input-file=$1.possorted.snap.nsrt.bam  \
	--output-snap=$1.possorted.snap  \
	--genome-name=hg38  \
	--genome-size=hg38.chrom.sizes  \
	--min-mapq=30  \
	--min-flen=50  \
	--max-flen=1000  \
	--keep-chrm=TRUE  \
	--keep-single=FALSE  \
	--keep-secondary=False  \
	--overwrite=True  \
	--max-num=20000  \
	--min-cov=500  \
	--verbose=True

snaptools snap-add-bmat \
    --snap-file=$1.possorted.snap \
    --bin-size-list 1000 5000 10000 \
    --verbose=True