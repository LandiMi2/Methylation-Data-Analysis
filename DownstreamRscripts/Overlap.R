#Per context. Modify script using sed command for other context

library(data.table)
library(GenomicRanges)
setwd("~/Documents/PhD-Bioinformatics-EpiCass/Analysis/Arab-virus-analysis/NewAnalysisNov2025/BismarkAnalysis/dmrs/")

dmr <- fread("mock_cmv/mock_cmv_cpg.0.10.dmrs",
             col.names = c("chr","start","end","qvalue","methDiff","CpN","MWU","2DKS",
                           "mock","camv"))

##get hypo and hyper
dmr.hyper <- dmr[methDiff>0]
dmr.hypo <- dmr[methDiff<0]

## read gene annotations
genes <- fread("../../glob/genes.bed", col.names = c("chr","start","end","id","blank","strand"))


#promoters
promoter <- genes[, .(chr,
                      start = ifelse(strand == "+", pmax(0, start - 2000), end + 1),
                      end = ifelse(strand == "+", start - 1, end + 2000), 
                      strand,id )]

### genomics ranges for DMRs  - hyper and hypo
dmr.gr.dmr.hyper <- GRanges(
  seqnames = dmr.hyper$chr,
  ranges = IRanges(start = dmr.hyper$start, 
                   end = dmr.hyper$end),
  diff.Methy = dmr.hyper$methDiff,
  nCGs = dmr.hyper$CpN
)

dmr.gr.dmr.hypo <- GRanges(
  seqnames = dmr.hypo$chr,
  ranges = IRanges(start = dmr.hypo$start, 
                   end = dmr.hypo$end),
  diff.Methy = dmr.hypo$methDiff,
  nCGs = dmr.hypo$CpN
)


######  get overlapping genes to CpG hyper and hypo
#gff to Granges
genes.gr <- GRanges(
  seqnames = genes$chr,
  ranges = IRanges(start = genes$start, 
                   end = genes$end),
  geneId = genes$id 
)

#get overlap for CpG for call cultivars
overlapsGeneDMRs.dmr.hyper <- findOverlaps(dmr.gr.dmr.hyper, genes.gr)
overlapsGeneDMRs.dmr.hypo <- findOverlaps(dmr.gr.dmr.hypo, genes.gr)


#extract the overlapping regions
overlap.dmr.dmr.hyper <-  dmr.gr.dmr.hyper[queryHits(overlapsGeneDMRs.dmr.hyper )]
overlap.genes.dmr.hyper <-  genes.gr[subjectHits(overlapsGeneDMRs.dmr.hyper )]

overlap.dmr.dmr.hypo <-  dmr.gr.dmr.hypo[queryHits(overlapsGeneDMRs.dmr.hypo)]
overlap.genes.dmr.hypo <-  genes.gr[subjectHits(overlapsGeneDMRs.dmr.hypo)]

##create an overlap data frame for all cultivars for all cultivars 
overlap.DMRgenes.dmr.hyper <- data.frame(
  DMRChr = seqnames(overlap.dmr.dmr.hyper),
  DMRsStart = start(overlap.dmr.dmr.hyper ),
  DMRsEnd = end(overlap.dmr.dmr.hyper ),
  DMRsCGCount = overlap.dmr.dmr.hyper$nCGs, 
  DMRsDiffMethy = overlap.dmr.dmr.hyper$diff.Methy,
  AnnoChr = seqnames(overlap.genes.dmr.hyper),
  AnnoStart = start(overlap.genes.dmr.hyper),
  AnnoEnd = end(overlap.genes.dmr.hyper),
  GeneID = overlap.genes.dmr.hyper$geneId)

setDT(overlap.DMRgenes.dmr.hyper)
overlap.DMRgenes.dmr.hyper

overlap.DMRgenes.dmr.hypo <- data.frame(
  DMRChr = seqnames(overlap.dmr.dmr.hypo),
  DMRsStart = start(overlap.dmr.dmr.hypo ),
  DMRsEnd = end(overlap.dmr.dmr.hypo),
  DMRsCGCount = overlap.dmr.dmr.hypo$nCGs, 
  DMRsDiffMethy = overlap.dmr.dmr.hypo$diff.Methy,
  AnnoChr = seqnames(overlap.genes.dmr.hypo),
  AnnoStart = start(overlap.genes.dmr.hypo),
  AnnoEnd = end(overlap.genes.dmr.hypo),
  GeneID = overlap.genes.dmr.hypo$geneId)

setDT(overlap.DMRgenes.dmr.hypo)
overlap.DMRgenes.dmr.hypo



###promoters ###############################
promoter.gr <- GRanges(
  seqnames = promoter$chr,
  ranges = IRanges(start = promoter$start, 
                   end = promoter$end),
  promterId = promoter$id 
)

#overlap for CpGs
overlapsPromoterDMRs.dmr.hyper <- findOverlaps(dmr.gr.dmr.hyper, promoter.gr)
overlapsPromoterDMRs.dmr.hypo <- findOverlaps(dmr.gr.dmr.hypo, promoter.gr)



#extract the overlapping regions
overlap.dmr.promoter.dmr.hyper <- dmr.gr.dmr.hyper[queryHits(overlapsPromoterDMRs.dmr.hyper)]
overlap.promoter.dmr.hyper <- promoter.gr[subjectHits(overlapsPromoterDMRs.dmr.hyper)]

overlap.dmr.promoter.dmr.hypo <- dmr.gr.dmr.hypo[queryHits(overlapsPromoterDMRs.dmr.hypo)]
overlap.promoter.dmr.hypo <- promoter.gr[subjectHits(overlapsPromoterDMRs.dmr.hypo)]


##create an overlap data frame
overlap.DMRpromoter.dmr.hyper <- data.frame(
  DMRChr = seqnames(overlap.dmr.promoter.dmr.hyper),
  DMRsStart = start(overlap.dmr.promoter.dmr.hyper),
  DMRsEnd = end(overlap.dmr.promoter.dmr.hyper),
  DMRsCGCount = overlap.dmr.promoter.dmr.hyper$nCGs,
  DMRsDiffMethy = overlap.dmr.promoter.dmr.hyper$diff.Methy,
  AnnoChr = seqnames(overlap.promoter.dmr.hyper),
  AnnoStart = start(overlap.promoter.dmr.hyper),
  AnnoEnd = end(overlap.promoter.dmr.hyper),
  GeneID = overlap.promoter.dmr.hyper$promterId)

setDT(overlap.DMRpromoter.dmr.hyper)
overlap.DMRpromoter.dmr.hyper

overlap.DMRpromoter.dmr.hypo <- data.frame(
  DMRChr = seqnames(overlap.dmr.promoter.dmr.hypo),
  DMRsStart = start(overlap.dmr.promoter.dmr.hypo),
  DMRsEnd = end(overlap.dmr.promoter.dmr.hypo),
  DMRsCGCount = overlap.dmr.promoter.dmr.hypo$nCGs,
  DMRsDiffMethy = overlap.dmr.promoter.dmr.hypo$diff.Methy,
  AnnoChr = seqnames(overlap.promoter.dmr.hypo),
  AnnoStart = start(overlap.promoter.dmr.hypo),
  AnnoEnd = end(overlap.promoter.dmr.hypo),
  GeneID = overlap.promoter.dmr.hypo$promterId)

setDT(overlap.DMRpromoter.dmr.hypo)
overlap.DMRpromoter.dmr.hypo

overlap.DMRgenes.dmr.hyper$context <- "CpG"
overlap.DMRgenes.dmr.hypo$context  <- "CpG"
overlap.DMRpromoter.dmr.hyper$context  <- "CpG"
overlap.DMRpromoter.dmr.hypo$context  <- "CpG"


#### write files
write.table(overlap.DMRgenes.dmr.hyper,"mock_cmv/mock_cmvDMRsGenes.CpG.Hyper.txt",sep="\t", row.names = F, quote = F)
write.table(overlap.DMRgenes.dmr.hypo,"mock_cmv/mock_cmvDMRsGenes.CpG.Hypo.txt",sep="\t", row.names = F, quote = F)

#promoter
write.table(overlap.DMRpromoter.dmr.hyper,"mock_cmv/mock_cmvDMRsPromoter.CpG.Hyper.txt",sep="\t", row.names = F, quote = F)
write.table(overlap.DMRpromoter.dmr.hypo,"mock_cmv/mock_cmvDMRsPromoter.CpG.Hypo.txt",sep="\t", row.names = F, quote = F)



