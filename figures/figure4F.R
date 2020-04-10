library(ggpubr)
library(BuenColors)
data("ToothGrowth")
ggbarplot(ToothGrowth, x = "dose", y = "len", 
          add = c("mean_se", "jitter"))
#vcam1
vcam1 <- read.csv("fig4g_data/vcam1.csv")
fig4g_vcam1 <- ggbarplot(vcam1, x = "TNFa.treatment", y = "Expression", 
               add = c("mean_sd", "jitter"),color = "TNFa.treatment", palette = c("#00AFBB", "#bb8d00", "#bb2600"),order = c("Control","24 h", "48 h"),
               position = position_dodge(0.8),ylim = c(0, 3))+scale_y_continuous(expand = c(0, 0))+stat_compare_means(comparisons = list(c("Control","24 h"),c("Control","48 h")), hide.ns = T,method = "t.test", label = "p.signif") #270x380

#tpm1
tpm1 <- read.csv("fig4g_data/tpm1.csv")
fig4g_tpm1  <- ggbarplot(tpm1, x = "TNFa.treatment", y = "Expression", 
          add = c("mean_sd", "jitter"),color = "TNFa.treatment", palette = c("#00AFBB", "#bb8d00", "#bb2600"),order = c("Control","24 h", "48 h"),
          position = position_dodge(0.8),ylim = c(0, 1.5))+scale_y_continuous(expand = c(0, 0))+stat_compare_means(comparisons = list(c("Control","48 h")),method = "t.test", label = "p.signif") #270x380

#vcam1
slc5a12 <- read.csv("fig4g_data/slc5a12.csv")
fig4g_slc5a12 <- ggbarplot(slc5a12, x = "TNFa.treatment", y = "Expression", 
                         add = c("mean_sd", "jitter"),color = "TNFa.treatment", palette = c("#00AFBB", "#bb8d00", "#bb2600"),order = c("Control","24 h", "48 h"),
                         position = position_dodge(0.8),ylim = c(0, 1.5))+scale_y_continuous(expand = c(0, 0))+stat_compare_means(comparisons = list(c("Control","48 h")),method = "t.test", label = "p.signif") #270x380

#tpm1
slc4a4 <- read.csv("fig4g_data/slc4a4.csv")
fig4g_slc4a4  <- ggbarplot(slc4a4, x = "TNFa.treatment", y = "Expression", 
                         add = c("mean_sd", "jitter"),color = "TNFa.treatment", palette = c("#00AFBB", "#bb8d00", "#bb2600"),order = c("Control","24 h", "48 h"),
                         position = position_dodge(0.8),ylim = c(0, 1.5))+scale_y_continuous(expand = c(0, 0))+stat_compare_means(comparisons = list(c("Control","24 h"),c("Control","48 h")),method = "t.test", label = "p.signif") #270x380