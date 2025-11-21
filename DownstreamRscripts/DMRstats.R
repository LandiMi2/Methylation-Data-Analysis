
setwd("~/Documents/PhD-Bioinformatics-EpiCass/Analysis/Arab-virus-analysis/NewAnalysisNov2025/BismarkAnalysis/dmrs/")

library(data.table)
library(ggplot2)

cpg <- fread("mock_cmv/mock_cmv_cpg.0.10.dmrs")
chg <- fread("mock_cmv/mock_cmv_chg.0.10.dmrs")
chh <- fread("mock_cmv/mock_cmv_chh.0.10.dmrs")

#get relevant colum 
cpg <- cpg[,.(V1,V2,V3,V5,V6,V9,V10)]
chg <- chg[,.(V1,V2,V3,V5,V6,V9,V10)]
chh <- chh[,.(V1,V2,V3,V5,V6,V9,V10)]

#column
colnames(cpg) <- c("chr","start","end","methDiff","sites","mock","camv")

#add hypo and hyper methylation 
cpg[, MethState := ifelse(methDiff > 0 , "hyper",
                          ifelse(methDiff < 0 , "hypo", 0))]

#chg
colnames(chg) <-  c("chr","start","end","methDiff","sites","mock","camv")
#add hypo and hyper methylation 
chg[, MethState := ifelse(methDiff > 0 , "hyper",
                          ifelse(methDiff < 0 , "hypo", 0))]
#
colnames(chh) <-  c("chr","start","end","methDiff","sites","mock","camv")
#add hypo and hyper methylation 
chh[, MethState := ifelse(methDiff > 0 , "hyper",
                          ifelse(methDiff < 0 , "hypo", 0))]



### plot the number of hyper and hypo methylation 
#create a context column 
cpg[, context := "CpG"]
chg[, context := "CHG"]
chh[, context := "CHH"]

#merge 
combined <- rbind(cpg, chg,chh)
#frequencies by context and class
freq.dat <- combined[, .N, by = .(context, MethState)]
freq.dat[, total := sum(N), by = context]
freq.dat[, percent := round((N / total) * 100)]
freq.dat[, label := paste0(N, " (", percent, "%)")]
#order
freq.dat[, context := factor(context, levels = c("CpG", "CHG", "CHH"))]

#
#plot 
ggplot(freq.dat, aes(x = MethState, y = N, fill = MethState)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = label), vjust = -0.5, size = 4) + 
  facet_wrap(~context, scales = "free_y") +
  labs(
    y = "DMRs",
    fill = "Methylation State",
    title = "Mock vs CMV"
  ) +
  theme_classic() +
  scale_fill_manual(values = c("hyper" = "#ff7f0e", "hypo" = "#1f77b4"))




