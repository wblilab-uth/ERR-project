# supplementary notes Fig 6
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

# binomial test -----------------------------------------------------------
# p = promoter mutation probability (success)
# N = total donors
# k =  dornor affected
# pbinom(q, size, prob, lower.tail = TRUE, log.p = FALSE)
# q = donor with the mutation
# size = donor of the cancer type
# prob = all mut / (gene * total donors) or mean(affected%)
# 37552 hg19.refGene.tss.uniq.srt.bed
# totalgene0 <- 27502
genepromoters <- fread("annotations/hg19.refGene.tss.count", header = F, sep = "\t")
cancerpmut.count <- as.data.table(table(unique(muts[,c(1,3,5)])[,1]))
mutinfo.2[, totalpMut := cancerpmut.count[match(mutinfo.2$Cancer, V1)]$N]
# mutinfo.2[, pmutexpect := totalpMut / 27502 / DonorCount]
mutinfo.2[, pmutexpect := totalpMut / 37552 / DonorCount]
mutinfo.1[, GenePromoters := genepromoters[match(mutinfo.1$Gene, V1)]$V2]
mutinfo.1[, pmutexpect := mutinfo.2[match(mutinfo.1$Cancer, Cancer)]$pmutexpect]
mutinfo.1[, hotspot.p := pbinom(N, size = CancerDonors, prob = pmutexpect * GenePromoters, lower.tail = F)]
# mutinfo.1[is.na(hotspot.p)]
qvalues <- p.adjust(mutinfo.1$hotspot.p, method = "fdr")
mutinfo.1 <- data.table(mutinfo.1, FDRs = qvalues)
sigonconeighbor.mut <- unique(mutinfo.1[FDRs < 1e-2 & Gene %in% selectedgenes$neighboringGene]$Gene)
sigoncogene.mut <- unique(mutinfo.1[FDRs < 1e-2 & Gene %in% selectedgenes$oncogene]$Gene)

cancerpdel.count <- as.data.table(table(unique(dels[,c(13,14,16)])[,1]))
delinfo.2[, totalpDel := cancerpdel.count[match(delinfo.2$V1, V1)]$N]
# mutinfo.2[, pmutexpect := totalpMut / 27502 / DonorCount]
delinfo.2[, pdelexpect := totalpDel / 37552 / V2]
delinfo.1[, GenePromoters := genepromoters[match(delinfo.1$Gene, V1)]$V2]
delinfo.1[, pdelexpect := delinfo.2[match(delinfo.1$Cancer, V1)]$pdelexpect]
delinfo.1[, hotspot.p := pbinom(N, size = CancerDonors, prob = pdelexpect * GenePromoters, lower.tail = F)]
# delinfo.1[is.na(hotspot.p)]
qvalues <- p.adjust(delinfo.1$hotspot.p, method = "fdr")
delinfo.1 <- data.table(delinfo.1, FDRs = qvalues)
# delinfo.1[hotspot.p < 1e-3]
# delinfo.1[FDRs < 1e-2][is.na(FDRs)]
# delinfo.1[Gene == "TP53"]
# delinfo.1[Gene == "TERT" & FDRs < 1e-2]
sigonconeighbor.del <- unique(delinfo.1[FDRs < 1e-2 & Gene %in% selectedgenes$neighboringGene]$Gene)
sigoncogene.del <- unique(delinfo.1[FDRs < 1e-2 & Gene %in% selectedgenes$oncogene]$Gene)

# show cancers with deletion < 1e6  -------------------------------------
selectedgenes.muts <- unique(muts[Gene %in% selectedgenes$neighboringGene | Gene %in% selectedgenes$oncogene, c(1,4)])
selectedgenes.dels <- unique(dels[DelLength < 1e6][Gene %in% selectedgenes$neighboringGene | Gene %in% selectedgenes$oncogene, c(13,16)])
selectedgenes.muts[, ID := paste(Cancer, Gene)]
selectedgenes.dels[, ID := paste(Cancer, Gene)]
mutinfo.1[, ID := paste(Cancer, Gene)]
delinfo.1[, ID := paste(Cancer, Gene)]
selectedgenes.muts[, SamplePerc := mutinfo.1[match(selectedgenes.muts$ID, ID)]$SamplePerc]
selectedgenes.dels[, SamplePerc := delinfo.1[match(selectedgenes.dels$ID, ID)]$SamplePerc]
selectedgenes.muts[, DonorCount := mutinfo.1[match(selectedgenes.muts$ID, ID)]$N]
selectedgenes.dels[, DonorCount := delinfo.1[match(selectedgenes.dels$ID, ID)]$N]
selectedgenes.muts[, CancerDonors := mutinfo.1[match(selectedgenes.muts$ID, ID)]$CancerDonors]
selectedgenes.dels[, CancerDonors := delinfo.1[match(selectedgenes.dels$ID, ID)]$CancerDonors]
selectedgenes.muts[, TotalDonorCount := sum(DonorCount), by = "Gene"]
selectedgenes.dels[, TotalDonorCount := sum(DonorCount), by = "Gene"]

selectedgenes.comb <- data.table(rbindlist(list(selectedgenes.muts, selectedgenes.dels)),
                                 Group = c(rep("mutation",nrow(selectedgenes.muts)),
                                           rep("deletion", nrow(selectedgenes.dels))))
selectedgenes.comb[, oncogene := selectedgenes[match(selectedgenes.comb$Gene, neighboringGene)]$oncogene]
selectedgenes.comb[, Labs := ifelse(is.na(oncogene), paste(Group, ":", Gene), paste(Group, ":", Gene, "->", oncogene))]
selectedgenes.comb[, Class := ifelse(is.na(oncogene), Gene, oncogene)]
selectedgenes.comb[, Class1 := ifelse(is.na(oncogene), "Oncogene", "Neighbor")]
setorder(selectedgenes.comb, Class, Labs)
selectedgenes.comb$Labs <- factor(selectedgenes.comb$Labs, levels = unique(selectedgenes.comb$Labs))
selectedgenes.comb$Cancer <- factor(selectedgenes.comb$Cancer, levels = cancersites$V1)
selectedgenes.comb$Class1 <- factor(selectedgenes.comb$Class1, levels = c("Neighbor", "Oncogene"))

p2 <- ggscatter(data = selectedgenes.comb[(Group == "mutation" & Gene %in% sigonconeighbor.mut) | 
                                              (Group == "deletion" & Gene %in% sigonconeighbor.del)],
                x = "Cancer", y = "Labs",
                size = "SamplePerc", 
                color = "Group", fill = "Group", palette = "aaas",
                # shape = "Group",
                xlab = "", ylab = "",
                title = "top oncoNeighbor genes (promoter: -1~1K, deletion: < 10Kb, FDR < 0.001 in >= 1 Cancer)") +
    scale_size(range = c(0.5, 5)) + theme_bw()
ggsave(filename = "ggscatter.toponcogene.mut.del.32.pdf", useDingbats = F, dpi = 300, 
       width = 15, height = 20, ggpar(p2, x.text.angle = 45))
