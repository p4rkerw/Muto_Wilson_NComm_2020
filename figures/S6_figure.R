library(ggpubr)
library(BuenColors)
#atac_overlap_ccans_gh.csv is extracted from atac_overlap_ccans_geneHancer.xlsx
ccan_gh <- read.csv("fig_csvdata/atac_overlap_ccans_gh.csv")
figs6 <- ggboxplot(ccan_gh, x = "threshold", y = "overlapping",
          add = "jitter",color = "threshold",ylim = c(0., 0.9))   