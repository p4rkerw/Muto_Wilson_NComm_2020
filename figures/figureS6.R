library(ggpubr)
library(BuenColors)

ccan_gh <- read.csv("fig_csvdata/atac_overlap_ccans_gh.csv")
ggboxplot(ccan_gh, x = "threshold", y = "overlapping",
          add = "jitter",color = "threshold",ylim = c(0., 0.9))     