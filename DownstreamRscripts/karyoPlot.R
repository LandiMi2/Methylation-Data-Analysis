### this R script plots per context, in this case CpG. Modify for other context. 
#input is from this command - multiBigwigSummary bins -b *.bw --binSize 100000 -out meth100kbCpG.npz --outRawCounts meth100kbCpG.tab 


library(data.table)
library(karyoploteR)
library(GenomicRanges)
setwd("~/Documents/PhD-Bioinformatics-EpiCass/Analysis/Arab-virus-analysis/NewAnalysisNov2025/BismarkAnalysis/glob/")

data <- fread("meth100kbCpG.tab")

gene <- fread("../../glob/genes.bed", col.names = c("chr","start","end","id","blank","strand"))

te <- fread("../../glob/TAIR10_GFF3_genes_transposons.gff")
te <- te[,.(V1,V3,V4,V5,V9)]
colnames(te) <- c("chr","anno","start","end","id")
te <- te[anno =="transposable_element"]
#make Chr column small c
te[, chr := tolower(chr)]

#genes
gene.gr <- GRanges(seqnames = gene$chr,
                   ranges = IRanges(start = gene$start, end = gene$end))

#te 
te.gr <- GRanges(seqnames = te$chr,
                 ranges = IRanges(start = te$start, end = te$end))



##for all 
colnames(data) <- c("chr","start","end","CMV1","CMV2","CMV3","CaMV1",
                    "CaMV2","CaMV3", "Mock1","Mock2","Mock3")

data[start == 0, start := 1]

## remove rep 1
data[, c("CMV1", "CaMV1", "Mock1") := NULL]


bins.gr <- GRanges(seqnames = data$chr,
                   ranges = IRanges(start = data$start, end = data$end))


#chromoses
chom <- fread("../../glob/chrom.sizes", header = F)
colnames(chom) <- c("chr", "end")
chom[, start := 1]
setcolorder(chom, c("chr", "start", "end"))
custom.genome <- toGRanges(chom)


samples <- c("CMV2","CMV3",
             "CaMV2","CaMV3",
             "Mock2","Mock3")

#color
colors <- c( "#ff2a0e","#ff2a0e",
             "#ff9f4e","#ff9f4e",
             "#1f77b4","#1f77b4")



pdf("CpG.pdf", width = 16, height = 14)

ordered.chrs <- paste0("chr", 1:5)

#plot parameters 
plot.params <- getDefaultPlotParams(plot.type = 1)
plot.params$data1outmargin <- 20  
plot.params$data1inmargin <- 10   
plot.params$ideogramheight <- 2  

kp <- plotKaryotype(
  genome = custom.genome,
  chromosomes = ordered.chrs,
  plot.type = 1,
  plot.params = plot.params
)

kpAddBaseNumbers(kp, tick.dist = 1e6, cex = 0.3)

## genes 
kpPlotDensity(kp,
              data = gene.gr,
              window.size = 100000, 
              col = "black",
              r0 = 0, r1 = 0.08)  

# tes
kpPlotDensity(kp,
              data = te.gr,
              window.size = 100000,
              col = "grey",  
              r0 = 0, r1 = 0.08)


####
for (i in seq_along(samples)) {
  meth.gr <- GRanges(seqnames = data$chr,
                     ranges = IRanges(start = data$start, end = data$end),
                     score = data[[samples[i]]])
  
  kpLines(kp,
          data = meth.gr,
          y = data[[samples[i]]] / 100,
          col = colors[i],
          lwd = 1,
          r0 = 0.1, r1 = 1)
}

kpAxis(kp, r0 = 0.1, r1 = 1, ymin = 0, ymax = 1, numticks = 2)

legend("topright",
       legend = c(samples, "Genes", "TEs"),
       col    = c(colors, "black", "grey"),
       lwd    = c(rep(2, length(samples)), 3, 3),
       border = NA,
       cex    = 0.8,
       bty    = "n",
       title  = "CpG")

dev.off()



