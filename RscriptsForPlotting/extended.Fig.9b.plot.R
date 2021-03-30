# screen -r errplot
# setwd("/data/shaojf/project.ERR/combined.gene.level")
library(data.table)
library(ggpubr)

`%notin%` <- Negate(`%in%`)
# read in general data ----------------------------------------------------------
cancersites <- fread("ICGC/ICGC.cancers", header = F)
setorder(cancersites, V3, V1)

geneannos <- fread("annotations/gencode.v19.annotation.gene.probemap", header = T, sep = "\t")
onconeighbors <- fread("Cosmic/allCosmic.200k.neigbors.gene.txt", header = F)
onconeighbors <- onconeighbors[V2 %notin% geneannos[chromEnd - chromStart < 1000]$gene][grep("MIR",V2, invert = T)]
onconeighbors <- unique(onconeighbors[V1 %notin% c("BIVM-ERCC5", "BUB1B-PAK6", "SEPT5-GP1BB", "HOXA10-HOXA9")])
colnames(onconeighbors) <- c("oncogene", "neighboringGene")
onconeighbors.tss <- fread("Cosmic/allCosmic.200k.neigbors.tss.txt", header = F)
oncogenes.tss <- fread("Cosmic/allCosmic.tss.txt", header = F)
cosmics <- fread("Cosmic/Census_allMon.Feb.10.17_05_43.2020.sim.tsv", header = T)
cosmics[grep(",", `Role in Cancer`), Role := "Multiple"]
cosmics[`Role in Cancer` == "oncogene", Role := "oncogene"]
cosmics[`Role in Cancer` == "TSG", Role := "TSG"]
cosmics[is.na(Role), Role := "Other"]

selectedgenes <- onconeighbors[oncogene %in% cosmics[grep("oncogene", `Role in Cancer`)]$`Gene Symbol`]
alloncogenes <- cosmics[grep("oncogene", `Role in Cancer`)]$`Gene Symbol`

muts <- fread("ICGC/ICGC.release28.WGS.refGene.tss.promoter.1k.mut.gene.dis.tss", header = T)
dels <- fread("ICGC/del.win1k.promoter.tsv", header = T)
colnames(dels)[13] <- "Cancer"

mutinfo.2 <- fread("ICGC/ICGC.release28.WGS.mut.donor.count.2", header = T)
mutinfo.1 <- as.data.table(table(unique(muts[,c(1,3:4)])[,c(1,3)]))
mutinfo.1 <- mutinfo.1[N > 0]
mutinfo.1[, CancerDonors := mutinfo.2[match(mutinfo.1$Cancer, Cancer)]$DonorCount]
mutinfo.1[, SamplePerc := N/CancerDonors * 100]

delinfo.2 <- fread("ICGC/ICGC.release28.WGS.del.donor.count.2", header = F)
delinfo.1 <- as.data.table(table(unique(dels[,c(13:14,16)])[,c(1,3)]))
delinfo.1 <- delinfo.1[N > 0]
delinfo.1[, CancerDonors := delinfo.2[match(delinfo.1$Cancer, V1)]$V2]
delinfo.1[, SamplePerc := N/CancerDonors * 100]

# mutinfo.1[, GenePromoters := genepromoters[match(mutinfo.1$Gene, V1)]$V2]
# delinfo.1[, GenePromoters := genepromoters[match(delinfo.1$Gene, V1)]$V2]

muts <- muts[grep("^MIR", Gene, invert = T)]
muts.count <- as.data.table(table(unique(muts[,c(3,5)])[,2]))

### merge all cancer type
allgenemut <- as.data.table(table(unique(muts[Mutation %in% muts.count[N>1]$V1,c(3,4)])[,2]))
allgenemut.1 <- as.data.table(table(unique(muts[Mutation %in% muts.count[N>2]$V1,c(3,4)])[,2]))
allgenemut.2 <- as.data.table(table(unique(muts[Mutation %in% muts.count[N>3]$V1,c(3,4)])[,2]))
allmutcount <- data.table(Cancer = c("ICGC.mut", "ICGC.mut>1", "ICGC.mut>2", "ICGC.mut>3", "ICGC.del<100K", "ICGC.del<10K"), 
                          N = c(length(unique(mutinfo.1$Gene)), 
                                length(unique(mutinfo.1[Gene %in% allgenemut$V1]$Gene)),
                                length(unique(mutinfo.1[Gene %in% allgenemut.1$V1]$Gene)),
                                length(unique(mutinfo.1[Gene %in% allgenemut.2$V1]$Gene)),
                                length(unique(dels[DelLength < 1e5]$Gene)),
                                length(unique(dels[DelLength < 1e4]$Gene))),
                          oncogene = c(length(unique(mutinfo.1[Gene %in% alloncogenes]$Gene)), 
                                       length(unique(mutinfo.1[Gene %in% allgenemut$V1 & Gene %in% alloncogenes]$Gene)),
                                       length(unique(mutinfo.1[Gene %in% allgenemut.1$V1 & Gene %in% alloncogenes]$Gene)),
                                       length(unique(mutinfo.1[Gene %in% allgenemut.2$V1 & Gene %in% alloncogenes]$Gene)),
                                       length(unique(dels[DelLength < 1e5 & Gene %in% alloncogenes]$Gene)),
                                       length(unique(dels[DelLength < 1e4 & Gene %in% alloncogenes]$Gene))),
                          neighboringGene = c(length(unique(mutinfo.1[Gene %in% selectedgenes$neighboringGene]$Gene)), 
                                              length(unique(mutinfo.1[Gene %in% allgenemut$V1 & Gene %in% selectedgenes$neighboringGene]$Gene)),
                                              length(unique(mutinfo.1[Gene %in% allgenemut.1$V1 & Gene %in% selectedgenes$neighboringGene]$Gene)),
                                              length(unique(mutinfo.1[Gene %in% allgenemut.2$V1 & Gene %in% selectedgenes$neighboringGene]$Gene)),
                                              length(unique(dels[DelLength < 1e5 & Gene %in% selectedgenes$neighboringGene]$Gene)),
                                              length(unique(dels[DelLength < 1e4 & Gene %in% selectedgenes$neighboringGene]$Gene)))
)
# total refgene: 27502
# totol cosmic oncogene: 315
# total cosmic onconeighbors: 1693
allmutcount[, oncogene.hyper.p := phyper(oncogene, N, 27502 - N, 315 - 1, lower.tail = FALSE, log.p = FALSE)]
allmutcount[, neighboringGene.hyper.p := phyper(neighboringGene, N, 27502 - N, 1693 - 1, lower.tail = FALSE, log.p = FALSE)]
allmutcount[, AffectedGenePercentage := N / 27502]
allmutcount[, AffectedOncoGenePercentage := oncogene / 315]
allmutcount[, AffectedOncoNeighborGenePercentage := neighboringGene / 1693]
x <- data.table(Cancer = allmutcount$Cancer, Group = rep("AffectedGene", nrow(allmutcount)), Percetange = 100 * allmutcount$AffectedGenePercentage, pvalue = NA)
y <- data.table(Cancer = allmutcount$Cancer, Group = rep("AffectedOncoGene", nrow(allmutcount)), Percetange = 100 * allmutcount$AffectedOncoGenePercentage, pvalue = allmutcount$oncogene.hyper.p)
z <- data.table(Cancer = allmutcount$Cancer, Group = rep("AffectedOncoNeighborGene", nrow(allmutcount)), Percetange = 100 * allmutcount$AffectedOncoNeighborGenePercentage, pvalue = allmutcount$neighboringGene.hyper.p)
allmutcount.melt <- rbindlist(list(x, y, z))
allmutcount.melt[pvalue < 1e-2, Labs := "*"]
allmutcount.melt[pvalue < 1e-3, Labs := "**"]
allmutcount.melt$Group <- factor(allmutcount.melt$Group, levels = c("AffectedGene", "AffectedOncoGene", "AffectedOncoNeighborGene"))
# allmutcount.melt[!is.na(Labs)]

p1 <- ggbarplot(data = allmutcount.melt[Cancer %in% c("ICGC.mut", "ICGC.mut>1", "ICGC.del<100K")], 
                x = "Cancer", y = "Percetange",
                color = "Group", fill = "Group", palette = "npg",
                width = 0.6, position = position_dodge(0.6), 
                label = as.vector(allmutcount.melt[Cancer %in% c("ICGC.mut", "ICGC.mut>1", "ICGC.del<100K")]$Labs), 
                lab.pos = "out", 
                xlab = "", ylab = "affected gene%",
                title = "affected genes (promoter: -1000~1000)")
ggsave(filename = "ggbarplot.promoter.affected.allstats.1.pdf", 
       useDingbats = F, dpi = 300, 
       width = 6, height = 6, ggpar(p1, x.text.angle = 0))
