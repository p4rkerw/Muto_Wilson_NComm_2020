# this script will identify the proportion cell type specific deg that also contain dar


# identify the genes associated with a differentially accessible chromatin region
dar_file <- "analysis_control/dar.celltype.control.xlsx"
dar_idents <- getSheetNames(dar_file)
dar.df.list <- lapply(dar_idents, function(x) read.xlsx(dar_file, sheet = x))
names(dar.df.list) <- dar_idents

dar.genes.list <- sapply(dar_idents, function(x) {
  genes <- dar.df.list[[x]] %>%
    dplyr::select(gene)
  return(genes)
})
names(dar.genes.list) <- dar_idents

# calculate total number of dar and number of unique dar genes
dar.total.num <- sapply(dar_idents, function(x) length(dar.genes.list[[x]]))
dar.uniqueGene.num <- sapply(dar_idents, function(x) length(unique(dar.genes.list[[x]])))


# identify the genes that are differentially expressed
deg_file <- "analysis_control/deg.celltype.control.xlsx"
deg_idents <- getSheetNames(deg_file)
deg.df.list <- lapply(deg_idents, function(x) {
  df <- read.xlsx(deg_file, sheet = x)
  colnames(df)[1] <- "gene" # add a column name
  return(df)
  })
names(deg.df.list) <- deg_idents

deg.genes.list <- sapply(deg_idents, function(x) {
  genes <- deg.df.list[[x]] %>%
    dplyr::select(gene)
  return(genes)
})
names(deg.genes.list) <- deg_idents

# calculate total number of deg
deg.total.num <- sapply(deg_idents, function(x) length(deg.genes.list[[x]]))


# prepare the dar idents for comparison with deg
deg_compare_idents <- c("PT","PT","PT_KIM1","PEC","TAL","DCT1","DCT2","CNT","PC","ICA","ICB","PODO","ENDO","MES",
                        "FIB","LEUK")
dar_compare_idents <- c("PCT","PST","PT_KIM1","PEC","TAL","DCT","DCT","CNT","PC","ICA","ICB","PODO","ENDO","MES_FIB",
                        "MES_FIB","LEUK")

# calculate number of genes with a deg and a unique dar for each celltype defined by deg
# ie calculate number of deg/dar overlaps
overlap.unique.list <- sapply(seq(deg_compare_idents), function(x){
  deg_ident <- deg_compare_idents[x]
  dar_ident <- dar_compare_idents[x]
  
  deg <- deg.genes.list[[deg_ident]]
  dar.unique <- unique(dar.genes.list[[dar_ident]])
  overlap <- intersect(deg, dar.unique)
  return(overlap)
})

names(overlap.unique.list) <- lapply(seq(overlap.unique.list), function(x) {
  if(deg_compare_idents[x] != dar_compare_idents[x]) {
  return(paste0(deg_compare_idents[x],"_vs_", dar_compare_idents[x]))
  } else {
  return(deg_compare_idents[x])
  }
})

overlap.unique.num.list <- sapply(names(overlap.unique.list), function(x) length(overlap.unique.list[[x]]))


# calculate total overlap for all dar (not just unique)
overlap.total.dar.list <- sapply(seq(deg_compare_idents), function(x){
  deg_ident <- deg_compare_idents[x]
  dar_ident <- dar_compare_idents[x]
  
  deg <- deg.genes.list[[deg_ident]]
  dar <- dar.genes.list[[dar_ident]]
  dar_no_deg <- setdiff(dar, deg) # list of dar not associated with a gene
  dar_with_deg <- dar[!(dar %in% dar_no_deg)] # note that setdiff will eliminate duplicate values
  return(dar_with_deg)
})

names(overlap.total.dar.list) <- lapply(seq(overlap.total.dar.list), function(x) {
  if(deg_compare_idents[x] != dar_compare_idents[x]) {
    return(paste0(deg_compare_idents[x],"_vs_", dar_compare_idents[x]))
  } else {
    return(deg_compare_idents[x])
  }
})

overlap.total.dar.num.list <- sapply(names(overlap.total.dar.list), function(x) length(overlap.total.dar.list[[x]]))






