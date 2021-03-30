# screen -r errplot
# setwd("/data/shaojf/project.ERR/combined.gene.level")
# the neighbors are not limited to 200k
library(data.table)
library(ggpubr)

`%notin%` <- Negate(`%in%`)
selectedgenes <- fread("otherlist/cripsr.genes", header = T)
# read in general data ----------------------------------------------------------
cancersites <- fread("ICGC/ICGC.cancers", header = F)
setorder(cancersites, V3, V1)

# geneannos <- fread("annotations/gencode.v19.annotation.gene.probemap", header = T, sep = "\t")
# onconeighbors <- fread("Cosmic/allCosmic.200k.neigbors.gene.txt", header = F)
# onconeighbors <- onconeighbors[V2 %notin% geneannos[chromEnd - chromStart < 1000]$gene][grep("MIR",V2, invert = T)]
# onconeighbors <- unique(onconeighbors[V1 %notin% c("BIVM-ERCC5", "BUB1B-PAK6", "SEPT5-GP1BB", "HOXA10-HOXA9")])
# colnames(onconeighbors) <- c("oncogene", "neighboringGene")
# onconeighbors.tss <- fread("Cosmic/allCosmic.200k.neigbors.tss.txt", header = F)
# oncogenes.tss <- fread("Cosmic/allCosmic.tss.txt", header = F)
# cosmics <- fread("Cosmic/Census_allMon.Feb.10.17_05_43.2020.sim.tsv", header = T)
# cosmics[grep(",", `Role in Cancer`), Role := "Multiple"]
# cosmics[`Role in Cancer` == "oncogene", Role := "oncogene"]
# cosmics[`Role in Cancer` == "TSG", Role := "TSG"]
# cosmics[is.na(Role), Role := "Other"]

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

####
selectedgenes.muts <- unique(muts[Gene %in% selectedgenes$neighboringGene | Gene %in% selectedgenes$oncogene, c(1,4)])
selectedgenes.muts[, ID := paste(Cancer, Gene)]
mutinfo.1[, ID := paste(Cancer, Gene)]
selectedgenes.muts[, SamplePerc := mutinfo.1[match(selectedgenes.muts$ID, ID)]$SamplePerc]
selectedgenes.muts[, DonorCount := mutinfo.1[match(selectedgenes.muts$ID, ID)]$N]

selectedgenes.dels <- unique(dels[DelLength < 1e5][Gene %in% selectedgenes$neighboringGene | Gene %in% selectedgenes$oncogene, c(13,16)])
selectedgenes.dels[, ID := paste(Cancer, Gene)]
delinfo.1[, ID := paste(Cancer, Gene)]
selectedgenes.dels[, SamplePerc := delinfo.1[match(selectedgenes.dels$ID, ID)]$SamplePerc]
selectedgenes.dels[, DonorCount := delinfo.1[match(selectedgenes.dels$ID, ID)]$N]

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

p1 <- ggscatter(data = selectedgenes.comb[Class1 == "Neighbor"],
                x = "Cancer", y = "Labs",
                size = "SamplePerc",
                color = "Group", fill = "Group", palette = "aaas",
                # shape = "Group",
                xlab = "", ylab = "",
                title = "selected Neighbor genes (promoter mutation: -1000~1000, deletion length: < 100Kb & DonorCount >= 1)") +
    scale_size(range = c(1, 8)) + theme_bw()
ggsave(filename = "ggscatter.selected.genes.mut.del.36pairs.pdf", useDingbats = F, dpi = 300,
       width = 12, height = 8, ggpar(p1, x.text.angle = 45))

p1 <- ggscatter(data = selectedgenes.comb[Class1 == "Neighbor" & DonorCount > 1],
                x = "Cancer", y = "Labs",
                size = "SamplePerc", 
                color = "Group", fill = "Group", palette = "aaas",
                # shape = "Group",
                xlab = "", ylab = "",
                title = "selected Neighbor genes (promoter mutation: -1000~1000, deletion length: < 100Kb & DonorCount > 1)") +
    scale_size(range = c(1, 8)) + theme_bw()
ggsave(filename = "ggscatter.selected.genes.mut.del.36pairs.rec.pdf", useDingbats = F, dpi = 300, 
       width = 12, height = 8, ggpar(p1, x.text.angle = 45))