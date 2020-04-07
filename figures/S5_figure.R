# this script will compare the cicero CCAN for each cell type to see which overlap with interactions
# in the GeneHancer database within 50kb of the top DAR
library(data.table)
library(dplyr)
library(here)
library(openxlsx)

list.clusterID <- getSheetNames("analysis_control/atac_overlap_ccans_geneHancer.xlsx")
list.overlap <- lapply(list.clusterID, function(sheetName) {
  read.xlsx("analysis_control/atac_overlap_ccans_geneHancer.xlsx",
  sheet= sheetName)
  })

prop.overlap.df <-  lapply(seq(list.overlap), function(x) {
  df <- dplyr::select(list.overlap[[x]], "prop_cicero_in_gh_de")
}) %>%
  bind_cols() 
colnames(prop.overlap.df) <- list.clusterID
rownames(prop.overlap.df) <- list(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)

# plot the results by transposing the dataframe so coaccess threshold are the columns 
boxplot(t(prop.overlap.df),
        main="Cicero connections overlap with GeneHancer 'Double Elite'",
        xlab="Cicero Coaccess Threshold",
        ylab="Mean Proportion of Cicero Connections in GeneHancer",
        col="red")
