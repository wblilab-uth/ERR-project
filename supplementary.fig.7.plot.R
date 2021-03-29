library(data.table)
library(ggpubr)

`%notin%` <- Negate(`%in%`)

# read in general data ----------------------------------------------------------
geneannos <- fread("annotations/gencode.v19.annotation.gene.probemap", header = T, sep = "\t")
colnames(geneannos)[1] <- "id"
cancersites <- fread("ICGC/ICGC.cancers", header = F)
setorder(cancersites, V3, V1)
totalmuts <- fread("ICGC/allstats.simple_somatic_mutation.open.tsv", header = F, sep = "\t")

onconeighbors <- fread("Cosmic/allCosmic.200k.neigbors.gene.txt", header = F)
onconeighbors <- onconeighbors[V2 %notin% geneannos[chromEnd - chromStart < 1000]$gene][grep("MIR",V2, invert = T)]
onconeighbors <- onconeighbors[V1 %notin% c("BIVM-ERCC5", "BUB1B-PAK6", "SEPT5-GP1BB", "HOXA10-HOXA9")]
onconeighbors.tss <- fread("Cosmic/allCosmic.200k.neigbors.tss.txt", header = F)
oncogenes.tss <- fread("Cosmic/allCosmic.tss.txt", header = F)
cosmics <- fread("Cosmic/Census_allMon.Feb.10.17_05_43.2020.sim.tsv", header = T)
cosmics[grep(",", `Role in Cancer`), Role := "Multiple"]
cosmics[`Role in Cancer` == "oncogene", Role := "oncogene"]
cosmics[`Role in Cancer` == "TSG", Role := "TSG"]
cosmics[is.na(Role), Role := "Other"]


ref.ctcf <- fread("ICGC/ICGC.release28.simple.mutation.ref.fimo.txt", header = F, sep = "\t")
mut.ctcf <- fread("ICGC/ICGC.release28.simple.mutation.mut.fimo.txt", header = F, sep = "\t")
ref.ctcf[, ID := paste(V3,V4,V5)]
mut.ctcf[, ID := paste(V3,V4,V5)]
ref.ctcf[, motif := ifelse(V8 < 1e-4, "Strong","Weak")]
mut.ctcf[, motif := ifelse(V8 < 1e-4, "Strong","Weak")]
# all.ctcf <- merge.data.table(x = ref.ctcf[V4 < 20 & V5 > 21], 
#                              y = mut.ctcf[V4 < 20 & V5 > 21], 
#                              by = "ID", all = TRUE)
all.ctcf <- merge(x = ref.ctcf[V4 < 20 & V5 > 21], 
                  y = mut.ctcf[V4 < 20 & V5 > 21],
                  by = "ID", all = TRUE)
all.ctcf[, Class := paste(motif.x, "->", motif.y)]
all.ctcf[, Mutation := matrix(unlist(strsplit(ID, split = " ")), byrow = T, ncol = 3)[,1]]
# all.ctcf[motif.x == motif.y & V7.y > 3.5]
all.ctcf.sim <- all.ctcf[motif.x != motif.y | (is.na(motif.x) & V7.y > 3.5) | (is.na(motif.y) & V7.x > 3.5),c(1,11,22,24:25)]
all.ctcf.sim[, MotifSeq := paste(V10.x, "->", V10.y)]
all.ctcf.mutcount <- table(unique(all.ctcf.sim[,4:5])$Mutation)
conflict.ids <- names(all.ctcf.mutcount[all.ctcf.mutcount>1])
tmpx <- unique(all.ctcf.sim[,4:5])
tmpy <- data.table(Class = rep("Conflict", length(conflict.ids)), Mutation = conflict.ids)
all.ctcf.labs <- rbindlist(list(tmpx[Mutation %notin% conflict.ids], tmpy))

mut2tss.counts.1 <- fread("ICGC/ICGC.release28.WGS.refGene.tss.promoter.10k.count.11", header = T, sep = "\t")
mut.count <- fread("ICGC/ICGC.release28.WGS.mut.count.10k.combined.txt", header = T)
mut.count.d <- mut.count[Group == "Donor"]
mut.count.p <- mut.count[Group == "Project"]
mut2tss.counts.1[, CosmicGene := cosmics[match(mut2tss.counts.1$Gene,`Gene Symbol`)]$`Gene Symbol`]
mut2tss.counts.1[, ClosetoCosmicGene := onconeighbors[match(mut2tss.counts.1$Gene,V2)]$V1]
mut2tss.counts <- mut2tss.counts.1
mut2tss.counts[, CosmicGeneRole := cosmics[match(mut2tss.counts$CosmicGene, `Gene Symbol`)]$`Role in Cancer`]
mut2tss.counts[, ClosetoCosmicGeneRole := cosmics[match(mut2tss.counts$ClosetoCosmicGene, `Gene Symbol`)]$`Role in Cancer`]
mut2tss.counts[, DonorCount := mut.count.d[match(mut2tss.counts$Mutation,V1)]$N]
mut2tss.counts[, ProjectCount := mut.count.p[match(mut2tss.counts$Mutation,V1)]$N]
mut2tss.counts[, CTCFp := all.ctcf.labs[match(mut2tss.counts$Mutation, Mutation)]$Class]
mut2tss.counts[CTCFp %in% c("Strong -> Weak", "Strong -> NA", "Weak -> NA"), CTCF := "Loss"]
mut2tss.counts[CTCFp %in% c("Weak -> Strong", "NA -> Strong", "NA -> Weak"), CTCF := "Gain"]
mut2tss.counts[CTCFp %in% c("Strong -> Weak", "Strong -> NA"), CTCF1 := "Loss"]
mut2tss.counts[CTCFp %in% c("Weak -> Strong", "NA -> Strong"), CTCF1 := "Gain"]

mut2tss.counts[!is.na(CosmicGeneRole), Class := CosmicGeneRole]
mut2tss.counts[!is.na(ClosetoCosmicGeneRole), Class := paste("close to", ClosetoCosmicGeneRole)]
mut2tss.counts[is.na(Class), Class := "Other"]
mut2tss.counts[, Class := sub(Class, pattern = ", fusion", replacement = "")]
mut2tss.counts[Class == "close to " | Class == "close to fusion" | Class == "" | Class == "fusion", Class := "Other"]
# mut2tss.counts[!is.na(CosmicGeneRole) & Class != "Other", Labs := CosmicGene]
# mut2tss.counts[!is.na(ClosetoCosmicGeneRole) & Class != "Other", Labs :=  paste(Gene, "->", ClosetoCosmicGene)]
# change to ->
mut2tss.counts[!is.na(CosmicGeneRole), Labs := CosmicGene]
mut2tss.counts[!is.na(ClosetoCosmicGeneRole), Labs :=  paste(Gene, "->", ClosetoCosmicGene)]
# mut2tss.counts[!is.na(Labs)]
# mut2tss.counts[!is.na(Labs) & !is.na(CTCF)]
mut2tss.counts[, CTCFmotif := CTCF]
mut2tss.counts[is.na(CTCFmotif), CTCFmotif := "noChange"]
mut2tss.counts[, CTCFmotif1 := CTCF1]
mut2tss.counts[is.na(CTCFmotif1), CTCFmotif1 := "noStrongMotifChange"]
# mut2tss.counts[CTCF == "Conflict"]
mut2tss.counts[CTCFp == "Conflict", CTCFmotif := "Conflict"]
mut2tss.counts[CTCFp == "Conflict", CTCFmotif := "Conflict"]
mut2tss.counts$CTCFmotif <- factor(mut2tss.counts$CTCFmotif, levels = c("noChange", "Gain", "Loss"))
mut2tss.counts$CTCFmotif1 <- factor(mut2tss.counts$CTCFmotif1, levels = c("noStrongMotifChange", "Gain", "Loss"))

tmpdata1 <- mut2tss.counts[Distance < 1000 & Distance > -1000 & Gene %in% c(onconeighbors[V1=="TERT"]$V2, "TERT")]
tmpdata1[is.na(CTCF), CTCF := "noChange"]
tmpdata1[is.na(CTCF1), CTCF1 := "noStrongChange"]

# matched expression ----------------------------------------------------------
# icgc.expr <- fread("../ICGC.release28/targets.all.exp.sim.tsv", header = F, sep = "\t")
icgc.expr <- fread("ICGC/targets.all.exp_seq.sim.tsv", header = F, sep = "\t")
icgc.expr[grep("TERT|ENST00000310581|ENSG00000164362|NM_198253", V4), V4 := "TERT"]
icgc.expr[grep("CLPTM1L|ENST00000320895|ENSG00000049656|NM_030782", V4), V4 := "CLPTM1L"]
tmp.expr <- icgc.expr[V4 == "TERT"]
tmp.expr[, zscore := scale(V5), by = "V2"]
tmp.muts <- totalmuts[V1 %in% tmpdata1$Mutation]
tmp.expr[V3 %in% tmp.muts[V1 %in% tmpdata1[CTCFmotif == "Loss" & Gene == "TERT"]$Mutation]$V4, TERTp := "TERT.Loss"]
tmp.expr[V3 %in% tmp.muts[V1 %in% tmpdata1[CTCFmotif == "Gain" & Gene == "TERT"]$Mutation]$V4, TERTp := "TERT.Gain"]
tmp.expr[V3 %in% tmp.muts[V1 %in% tmpdata1[CTCFmotif == "Loss" & Gene == "CLPTM1L"]$Mutation]$V4, Neighborp := "CLPTM1L.Loss"]
tmp.expr[V3 %in% tmp.muts[V1 %in% tmpdata1[CTCFmotif == "Loss" & Gene == "LINC01511"]$Mutation]$V4, Neighborp := "LINC01511.Loss"]
tmp.expr[V3 %in% tmp.muts[V1 %in% tmpdata1[CTCFmotif == "noChange" & Gene == "TERT"]$Mutation]$V4, TERTp := "TERT.noChange"]
tmp.expr[V3 %in% tmp.muts[V1 %in% tmpdata1[CTCFmotif == "noChange" & Gene != "TERT"]$Mutation]$V4, Neighborp := "Neighbor.noChange"]
tmp.expr[V3 %notin% totalmuts$V4, MutStatus := "Unknown"]
tmp.expr[V3 %in% totalmuts$V4 & is.na(TERTp) & is.na(Neighborp), MutStatus := "WT WT"]
tmp.expr[V3 %in% tmp.muts[V1 %in% c("MU832963")]$V4, TERTp := "HotSpot-66.MU832963"]
tmp.expr[V3 %in% tmp.muts[V1 %in% c("MU830690")]$V4, TERTp := "HotSpot-88.MU830690"]
tmp.expr[is.na(MutStatus), MutStatus := ifelse(is.na(TERTp), Neighborp, TERTp)]
tmp.expr[, Class := MutStatus]
tmp.expr[MutStatus == "Neighbor.noChange" | MutStatus == "TERT.noChange", Class := "otherPromMut"]
tmp.expr[MutStatus %in% c("LINC01511.Loss", "CLPTM1L.Loss"), Class := "Neighbor.Loss"]

tmp.expr[V3 %in% tmp.muts[V1 %in% tmpdata1[Gene == "TERT"]$Mutation]$V4, TERTpm := "TERT.p.mt"]
tmp.expr[V3 %in% tmp.muts[V1 %in% tmpdata1[Gene == "CLPTM1L"]$Mutation]$V4, Neighborpm := "CLPTM1L.p.mt"]
tmp.expr[V3 %in% tmp.muts[V1 %in% tmpdata1[Gene != "TERT" & Gene != "CLPTM1L"]$Mutation]$V4, Neighborpm := "Other.p.mt"]
tmp.expr[V3 %notin% totalmuts$V4, MutStatus1 := "Unknown"]
tmp.expr[V3 %in% totalmuts$V4 & is.na(TERTpm) & is.na(Neighborpm), MutStatus1 := "WT WT"]
tmp.expr[V3 %in% tmp.muts[V1 %in% c("MU832963")]$V4, TERTpm := "HotSpot-66.MU832963"]
tmp.expr[V3 %in% tmp.muts[V1 %in% c("MU830690")]$V4, TERTpm := "HotSpot-88.MU830690"]
tmp.expr[is.na(MutStatus1), MutStatus1 := ifelse(is.na(TERTpm), Neighborpm, TERTpm)]
tmp.expr[, Class1 := MutStatus1]
tmp.expr$Class1 <- factor(tmp.expr$Class1, levels = c("CLPTM1L.p.mt", "Others", "HotSpot-66.MU832963", "HotSpot-88.MU830690"))

tmp.expr[, TotalTMM := mean(V5, trim = 0.2, na.rm = T), by = "V2"]
tmp.expr[, FCtoTotalTMM := V5/TotalTMM]
p2 <- ggboxplot(data = tmp.expr[V2 %in% c("PBCA-US","SKCM-US","LICA-FR")],
                x = "V2", y = "FCtoTotalTMM",
                ylim = c(0,5), palette = "npg", width = 0.6,
                label = "V3", font.label = list(color = "Class1"),
                repel = T, outlier.shape = NA,
                label.select = list(criteria = "`Class1` %in% c('CLPTM1L.p.mt')"),
                add = "jitter", add.params = list(fill = "Class1", color = "Class1", size = 1.5),
                xlab = "", ylab = "TERT", title = "RNA-Seq FoldChange to TotalTMM") + 
    font("legend.text", size = 3.2)
ggsave(filename = "ggboxplot.ICGC.release28.all.TERT.expr.4class.pdf", width = 4, height = 4, p2)



## manually
# 36/ (2096+69 - 36) : 130 / (4386 + 204-130) : 1470 / (68923 + 2561 - 1470)
a <- 130
b <- 4386 + 204-130
c <- 36
d <- 2096 + 69 - 36
e <- 1470
f <- 68923 + 2561 - 1470
# totalonco <- nrow(mut2tss.counts[Distance < 1000 & Distance > -1000 & DonorCount > 1 & Class %in% c("oncogene", "oncogene, TSG")])
# totalonp <- nrow(mut2tss.counts[Distance < 1000 & Distance > -1000 & DonorCount > 1 & Class %in% c("close to oncogene", "close to oncogene, TSG")])
# totalother <- nrow(mut2tss.counts[Distance < 1000 & Distance > -1000 & DonorCount > 1 & is.na(CosmicGene) & is.na(ClosetoCosmicGene)][grep("MIR",Gene, invert = T)])
# a <- nrow(mut2tss.counts[Distance < 1000 & Distance > -1000 & DonorCount > 1 & !is.na(CTCF) & Class %in% c("close to oncogene", "close to oncogene, TSG")][grep("MIR",Gene, invert = T)])
# b <- nrow(mut2tss.counts[Distance < 1000 & Distance > -1000 & DonorCount > 1 & is.na(CTCF) & Class %in% c("close to oncogene", "close to oncogene, TSG")][grep("MIR",Gene, invert = T)])
# c <- nrow(mut2tss.counts[Distance < 1000 & Distance > -1000 & DonorCount > 1 & !is.na(CTCF) & Class %in% c("oncogene", "oncogene, TSG")])
# d <- nrow(mut2tss.counts[Distance < 1000 & Distance > -1000 & DonorCount > 1 & is.na(CTCF) & Class %in% c("oncogene", "oncogene, TSG")])
# e <- nrow(mut2tss.counts[Distance < 1000 & Distance > -1000 & DonorCount > 1 & !is.na(CTCF) & is.na(CosmicGene) & is.na(ClosetoCosmicGene)][grep("MIR",Gene, invert = T)])
# f <- nrow(mut2tss.counts[Distance < 1000 & Distance > -1000 & DonorCount > 1 & is.na(CTCF) & is.na(CosmicGene) & is.na(ClosetoCosmicGene)][grep("MIR",Gene, invert = T)])

totalonco <- nrow(mut2tss.counts[Distance < 1000 & Distance > -1000 & DonorCount > 1 & Class %in% c("oncogene", "oncogene, TSG")][grep("MIR",Gene, invert = T)])
totalonp <- nrow(mut2tss.counts[Distance < 1000 & Distance > -1000 & DonorCount > 1 & Class %in% c("close to oncogene", "close to oncogene, TSG")][grep("MIR",Gene, invert = T)])
# totalother <- nrow(mut2tss.counts[Distance < 1000 & Distance > -1000 & DonorCount > 1 & is.na(CosmicGene) & is.na(ClosetoCosmicGene)][grep("MIR",Gene, invert = T)])
a <- nrow(mut2tss.counts[Distance < 1000 & Distance > -1000 & DonorCount > 1 & CTCFmotif == "Loss" & Class %in% c("close to oncogene", "close to oncogene, TSG")][grep("MIR",Gene, invert = T)])
b <- totalonp - a
c <- nrow(mut2tss.counts[Distance < 1000 & Distance > -1000 & DonorCount > 1 & CTCFmotif == "Loss" & Class %in% c("oncogene", "oncogene, TSG")])
d <- totalonco - c
e <- nrow(mut2tss.counts[Distance < 1000 & Distance > -1000 & DonorCount > 1 & CTCFmotif == "Loss"])
f <- nrow(mut2tss.counts[Distance < 1000 & Distance > -1000 & DonorCount > 1][grep("MIR",Gene, invert = T)]) - e 

# totalonco <- nrow(mut2tss.counts[Distance < 1000 & Distance > -1000 & DonorCount > 1 & Class %in% c("oncogene", "oncogene, TSG")])
# totalonp <- nrow(mut2tss.counts[Distance < 1000 & Distance > -1000 & DonorCount > 1 & Class %in% c("close to oncogene", "close to oncogene, TSG")])
# totalother <- nrow(mut2tss.counts[Distance < 1000 & Distance > -1000 & DonorCount > 1 & is.na(CosmicGene) & is.na(ClosetoCosmicGene)])
# a <- nrow(mut2tss.counts[Distance < 1000 & Distance > -1000 & DonorCount > 1 & !is.na(CTCF) & Class %in% c("close to oncogene", "close to oncogene, TSG")])
# b <- nrow(mut2tss.counts[Distance < 1000 & Distance > -1000 & DonorCount > 1 & is.na(CTCF) & Class %in% c("close to oncogene", "close to oncogene, TSG")])
# c <- nrow(mut2tss.counts[Distance < 1000 & Distance > -1000 & DonorCount > 1 & !is.na(CTCF) & Class %in% c("oncogene", "oncogene, TSG")])
# d <- nrow(mut2tss.counts[Distance < 1000 & Distance > -1000 & DonorCount > 1 & is.na(CTCF) & Class %in% c("oncogene", "oncogene, TSG")])
# e <- nrow(mut2tss.counts[Distance < 1000 & Distance > -1000 & DonorCount > 1 & !is.na(CTCF) & is.na(CosmicGene) & is.na(ClosetoCosmicGene)])
# f <- nrow(mut2tss.counts[Distance < 1000 & Distance > -1000 & DonorCount > 1 & is.na(CTCF) & is.na(CosmicGene) & is.na(ClosetoCosmicGene)])

# "all promoter, not only other"
fisher.p4 <- c(fisher.test(matrix(c(a,b,c,d), nrow = 2))$p.value, fisher.test(matrix(c(a,b,e,f), nrow = 2))$p.value)
promoter.ctcf4 <- data.table(Count = c(a,b,c,d,e,f), 
                             Ratio = c(a/(a+b),b/(a+b),c/(c+d),d/(c+d),e/(e+f),f/c(e+f)),
                             CTCF = rep(c("CTCF-loss", "CTCF-noloss"), 3), 
                             Gene = c(rep("OncoNeighbors",2), rep("OncoGenes",2), rep("Others",2)))
promoter.ctcf4$CTCF <- factor(promoter.ctcf4$CTCF, levels = c("CTCF-loss", "CTCF-noloss"))
promoter.ctcf4$Gene <- factor(promoter.ctcf4$Gene, levels = c("OncoGenes", "OncoNeighbors", "Others"))
stat.test4 <- data.table(group1 = c("OncoNeighbors", "OncoNeighbors"), 
                         group2 = c("OncoGenes", "Others"), 
                         p = round(fisher.p4, 3), y.position = c(1.05,1.05))
p5 <- ggbarplot(promoter.ctcf4, x = "Gene", y = "Ratio",
                fill = "CTCF", color = "CTCF", palette = "npg", width = 0.6,
                label = "Count", lab.vjust = 1.1, lab.size = 3,
                xlab = "", ylab = "Ratio",
                title = "-1~1Kb & DonorCount > 1") + stat_pvalue_manual(stat.test4, label = "p")
ggsave(filename = "ggbarplot.ICGC.release28.WGS.refGene.tss.promoter.1k.CTCF.5.pdf", 
       width = 4, height = 5, p5)
