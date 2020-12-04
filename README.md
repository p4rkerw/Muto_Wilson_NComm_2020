# Muto_Wilson_bioRxiv_2020
Integration of paired snRNAseq and snATACseq from 5 healthy adult kidney cortex samples

Welcome to our github repository!  
Here you will find analysis scripts for our manuscript deposited in bioRxiv     
https://www.biorxiv.org/content/10.1101/2020.06.14.151167v1  

Please contact the co-first authors or corresponding author with questions or comments and visit the Humphrey's lab website at www.humphreyslab.com  
  
Thank you,  
Parker and Yoshi

Sample analysis and processing workflow
1. Generate a custom pre-mRNA index for cellranger (snRNA_prep)  
Libraries were generated from a nuclear dissociation and require a custom pre-mRNA reference to count introns. We used refdata-cellranger-GRCh38-3.0.0 which can be downloaded from the 10X genomics website: https://support.10xgenomics.com/ . The gtf file is processed with cellranger_rna_mkref.sh to create the GRCh38-1.2.0_premrna reference.  

2. Count each of the five snRNA libraries with cellranger and GRCh38-1.2.0 (snRNA_prep)  
cellranger_rna_count.sh  

3. Aggregate the five snRNA libraries using the cellranger_rna_aggr.csv file (snRNA_prep)    
cellranger_rna_aggr.sh  

4. Process the aggregated snRNA library and exclude doublets identified with DoubletFinder. This script will generate a processed snRNA library file rnaAggr_control.rds (snRNA_prep)    
seurat_rna_process.sh  

5. Count each of the five snATAC libraries with cellranger-atac and refdata-cellranger-atac-GRCh38-1.2.0 (snATAC_prep)  
cellranger_atac_count.sh

6. Aggregate the five snATAC libraries using the cellranger_atac_aggr.csv file (snATAC_prep)      
cellranger_atac_aggr.sh  

7. Processs the aggregated snATAC library and integrate with the snRNA library to remove doublets with label transfer. This script will generate a processed snATAC library file atacAggr_control.rds (snATAC_prep)  
signac_atac_process.sh  

8. Identify cell-specific differentially expressed genes in the snRNA dataset (analysis)
find_rna_deg.R  

9. Identify cell-specific differentially accessible chromatin in the snATAC dataset (analysis)  
find_atac_dar.R  

10. Find differentially expressed genes with nearby differentially accessible chromatin regions (analysis)  
find_overlap_deg_dar.R  

11.

