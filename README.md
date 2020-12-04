# Muto_Wilson_bioRxiv_2020  
**Single cell transcriptional and chromatin accessibility profiling redefine cellular heterogeneity in the adult human kidney**    
*Yoshiharu Muto, *Parker C. Wilson, Haojia Wu, Sushrut S. Waikar, Benjamin D. Humphreys  
*These authors contributed equally  
doi: https://doi.org/10.1101/2020.06.14.151167  


Welcome to our github repository!  
Here you will find analysis scripts for our manuscript deposited in bioRxiv where we integrate paired snRNAseq and snATACseq from 5 healthy adult kidney cortex samples. Please contact the co-first authors or corresponding author with questions or comments.  
<br/>
![alt text](http://humphreyslab.com/wp-content/uploads/2015/12/favicon-H.jpg)  
Visit the Humphrey's lab website:   
www.humphreyslab.com  
<br/>
Check out our interactive datasets with Kidney Interactive mulTiomics (KIT):  
http://humphreyslab.com/SingleCell/
<br/><br/>
Find us on Twitter: 
<br/>
<br/>
  <a href="https://twitter.com/HumphreysLab?ref_src=twsrc%5Etfw" class="twitter-follow-button" data-show-count="false">Follow @HumphreysLab</a>
  
Thanks,  
Parker and Yoshi

**Sample analysis and processing workflow**  
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

8. Find cell-specific differentially expressed genes in the snRNA dataset (analysis)  
find_rna_deg.R  

9. Find cell-specific differentially accessible chromatin in the snATAC dataset (analysis)  
find_atac_dar.R

10. Compare snATAC peaks with dnase hypersensitive sites (analysis)  
find_atac_overlap_dnase.R  

11. Find differentially expressed genes with nearby differentially accessible chromatin regions (analysis)  
find_overlap_deg_dar.R  

12. Find cis-coaccessibility networks in the snATAC dataset with Cicero (analysis)  
find_atac_ccans.R  

13. Annotate cis-coaccessibility networks and create circos plots to visualize links (analysis)  
annotate_ccans_circos.R

14. Annotate cis-coaccessibility networks to find promoter-enhancer links (analysis)  
annotate_peaks_fantom.R

15. Annotate cis-coaccessibility networks to find overlap with the GeneHancer database (analysis)
get_atac_dar_genehancer_conns.R
find_atac_overlap_ccans_genehancer.R

16. Find transcription factor motif activity in the snATAC dataset with chromVAR (analysis)  
find_atac_chromVAR.R  

17. Correlate transcription factor motif activity in the snATAC dataset with gene expression in the snRNA dataset (analysis)  
corr_chromVar_TF_exp.R  

18. Correlate cicero gene activity in the snATAC dataset with gene expression in the snRNA dataset (analysis)  
corr_cicero_gene_exp.R  

19. Perform pseudotemporal ordering of the distal nephron with Monocle (analysis)  
pseudotime_distal_nephron.R  

**Deconvolution:**    
Each dataset is first prepared as an ExpressionSet and then deconvolved with Bisque using the snRNA library. When necessary, transcript counts are prepared with Salmon and GRCh38.
1. Fan et al Human Diabetic Nephropathy (PMID:31578193, GSE142025)  
salmon_count_fan.R  
find_bulk_degs_fan.R  
deconvolution_fan_BisqueRNA.R  

2. Liu et al Mouse Kidney IRI (PMID:24569379, GSE98622)  
find_bulk_degs_liu.R  
deconvolution_liu_BisqueRNA.R  

3. TCGA non-tumor Kidney  
find_bulk_degs_tcga.R  
deconvolution_tcga_BisqueRNA.R  

**Allele Specific Analysis:**    
These scripts can be run in publicly-available docker containers found at https://hub.docker.com/ username: p4rkerw  
Each script has an example command to run the corresponding docker container  

(Follow the steps in order) 
1. Genotype the snRNA or snATAC libraries using GATK (or obtain a vcf from another method)    
2. Annotate the genotyped vcf with GATK Funcotator to evaluate gnomAD MAF and variant context  
3. Process the variant annotation into a csv file  
4. Filter the genotyped vcf for variants that overlap coding transcripts and introns  
5. Apply the WASP pipeline to cellranger-aligned bam files and realign overlapping variants with STAR  
6. Get allele-specific counts with GATK ASEReadCounter  
7. Filter heterozygous SNV and perform allele-specific analysis with ASEP  

**Utility Scripts:**    
getContigLengths.py  
make_ucsc_tracks.R  

**Figures:**      
Code for generating figures in the manuscript




