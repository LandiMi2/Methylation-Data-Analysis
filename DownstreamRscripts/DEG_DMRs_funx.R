#Per context. Modify the script using the sed command for other context

library(data.table)

setwd("~/Documents/PhD-Bioinformatics-EpiCass/Analysis/Arab-virus-analysis/NewAnalysisNov2025/BismarkAnalysis/dmrs/mock_cmv/")

deg <- fread("DEG_cmv.csv")

cpg.hyper.gene <- fread("mock_cmvDMRsGenes.CpG.Hyper.txt")
cpg.hypo.gene <- fread("mock_cmvDMRsGenes.CpG.Hypo.txt")

cpg.hyper.prom <- fread("mock_cmvDMRsPromoter.CpG.Hyper.txt")
cpg.hypo.prom <- fread("mock_cmvDMRsPromoter.CpG.Hypo.txt")


#check overlappin DMR-genes
cpg.hyper.gene.rna <- merge(cpg.hyper.gene, 
                                deg[, .(Gene, Regulation)],
                                by.x = "GeneID",
                                by.y = "Gene",
                                all.x = F)

cpg.hypo.gene.rna <- merge(cpg.hypo.gene, 
                            deg[, .(Gene, Regulation)],
                            by.x = "GeneID",
                            by.y = "Gene",
                            all.x = F)


#promoter
cpg.hyper.prom.rna <- merge(cpg.hyper.prom, 
                            deg[, .(Gene, Regulation)],
                            by.x = "GeneID",
                            by.y = "Gene",
                            all.x = F)

cpg.hypo.prom.rna <- merge(cpg.hypo.prom, 
                           deg[, .(Gene, Regulation)],
                           by.x = "GeneID",
                           by.y = "Gene",
                           all.x = F)


write.table(cpg.hyper.gene.rna,"DEG_DMRs/CpG/cpg.hyper.gene.rna.txt",sep="\t", row.names = F, quote = F)
write.table(cpg.hypo.gene.rna,"DEG_DMRs/CpG/cpg.hypo.gene.rna.txt",sep="\t", row.names = F, quote = F)
write.table(cpg.hyper.prom.rna,"DEG_DMRs/CpG/cpg.hyper.prom.rna.txt",sep="\t", row.names = F, quote = F)
write.table(cpg.hypo.prom.rna ,"DEG_DMRs/CpG/cpg.hypo.prom.rna.txt",sep="\t", row.names = F, quote = F)

############# add functions
#load annotation file
anno <- fread("TAIR10_fun.txt", quote = "") #downloaded from https://www.arabidopsis.org/
anno <- anno[, .(Model_name, Short_description)]
colnames(anno) <- c("Gene","Function")
anno[, Gene := sub("\\..*", "", Gene)] #modify gene IDs


## add funtion
cpg.hyper.gene.rna <- fread("mock_cmv/DEG_DMRs/CpG/cpg.hyper.gene.rna.txt")
cpg.hypo.gene.rna <- fread("mock_cmv/DEG_DMRs/CpG/cpg.hypo.gene.rna.txt")
cpg.hyper.prom.rna <- fread("mock_cmv/DEG_DMRs/CpG/cpg.hyper.prom.rna.txt")
cpg.hypo.prom.rna <- fread("mock_cmv/DEG_DMRs/CpG/cpg.hypo.prom.rna.txt")


#find overapping
cpg.hyper.gene.rna <- merge(cpg.hyper.gene.rna, anno,
                            by.x = "GeneID",
                            by.y = "Gene",
                            all.x = T)

cpg.hypo.gene.rna <- merge(cpg.hypo.gene.rna, anno,
                            by.x = "GeneID",
                            by.y = "Gene",
                            all.x = T)

cpg.hyper.prom.rna <- merge(cpg.hyper.prom.rna, anno,
                            by.x = "GeneID",
                            by.y = "Gene",
                            all.x = T)

cpg.hypo.prom.rna <- merge(cpg.hypo.prom.rna, anno,
                           by.x = "GeneID",
                           by.y = "Gene",
                           all.x = T)


write.table(cpg.hyper.gene.rna,"mock_cmv/DEG_DMRs/CpG/funx/cpg.hyper.gene.rna.funx.txt",sep="\t", row.names = F, quote = F)
write.table(cpg.hypo.gene.rna,"mock_cmv/DEG_DMRs/CpG/funx/cpg.hypo.gene.rna.funx.txt",sep="\t", row.names = F, quote = F)
write.table(cpg.hyper.prom.rna,"mock_cmv/DEG_DMRs/CpG/funx/cpg.hyper.prom.rna.funx.txt",sep="\t", row.names = F, quote = F)
write.table(cpg.hypo.prom.rna ,"mock_cmv/DEG_DMRs/CpG/funx/cpg.hypo.prom.rna.funx.txt",sep="\t", row.names = F, quote = F)







