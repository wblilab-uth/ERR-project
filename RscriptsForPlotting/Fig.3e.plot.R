library(data.table)
library(ggpubr)

`%notin%` <- Negate(`%in%`)
# read in general data ----------------------------------------------------------
geneannos <- fread("annotations/gencode.v19.annotation.gene.probemap", header = T, sep = "\t")
onconeighbors <- fread("Cosmic/allCosmic.200k.neigbors.gene.txt", header = F)
onconeighbors <- onconeighbors[V2 %notin% geneannos[chromEnd - chromStart < 1000]$gene][grep("MIR",V2, invert = T)]
# onconeighbors <- unique(onconeighbors[V1 %notin% c("BIVM-ERCC5", "BUB1B-PAK6", "SEPT5-GP1BB", "HOXA10-HOXA9")])
colnames(onconeighbors) <- c("oncogene", "neighboringGene")
onconeighbors.tss <- fread("Cosmic/allCosmic.200k.neigbors.tss.txt", header = F)
oncogenes.tss <- fread("Cosmic/allCosmic.tss.txt", header = F)
cosmics <- fread("Cosmic/Census_allMon.Feb.10.17_05_43.2020.sim.tsv", header = T)
cosmics[grep(",", `Role in Cancer`), Role := "Multiple"]
cosmics[`Role in Cancer` == "oncogene", Role := "oncogene"]
cosmics[`Role in Cancer` == "TSG", Role := "TSG"]
cosmics[is.na(Role), Role := "Other"]

# specific locus & CTCF
ref.ctcf <- fread("ICGC/ICGC.release28.simple.mutation.ref.fimo.txt", header = F, sep = "\t")
mut.ctcf <- fread("ICGC/ICGC.release28.simple.mutation.mut.fimo.txt", header = F, sep = "\t")
ref.ctcf[, ID := paste(V3,V4,V5)]
mut.ctcf[, ID := paste(V3,V4,V5)]

ref.ctcf[, motif := ifelse(V8 < 1e-4, "Strong","Weak")]
mut.ctcf[, motif := ifelse(V8 < 1e-4, "Strong","Weak")]
all.ctcf <- merge.data.table(x = ref.ctcf[V4 < 20 & V5 > 21], 
                             y = mut.ctcf[V4 < 20 & V5 > 21], 
                             by = "ID", all = TRUE)
all.ctcf[, Class := paste(motif.x, "->", motif.y)]
all.ctcf[, Mutation := matrix(unlist(strsplit(ID, split = " ")), byrow = T, ncol = 3)[,1]]
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
mut2tss.counts.1[, ClosetoCosmicGene := onconeighbors[match(mut2tss.counts.1$Gene,neighboringGene)]$oncogene]
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
# change to ->
mut2tss.counts[!is.na(CosmicGeneRole), Labs := CosmicGene]
mut2tss.counts[!is.na(ClosetoCosmicGeneRole), Labs :=  paste(Gene, "->", ClosetoCosmicGene)]
mut2tss.counts[, CTCFmotif := CTCF]
mut2tss.counts[is.na(CTCFmotif), CTCFmotif := "noChange"]
mut2tss.counts[, CTCFmotif1 := CTCF1]
mut2tss.counts[is.na(CTCFmotif1), CTCFmotif1 := "noStrongMotifChange"]
mut2tss.counts[CTCFp == "Conflict", CTCFmotif := "Conflict"]
mut2tss.counts[CTCFp == "Conflict", CTCFmotif := "Conflict"]

mut2tss.counts$CTCFmotif <- factor(mut2tss.counts$CTCFmotif, levels = c("noChange", "Gain", "Loss"))
mut2tss.counts$CTCFmotif1 <- factor(mut2tss.counts$CTCFmotif1, levels = c("noStrongMotifChange", "Gain", "Loss"))
# mut2tss.counts[Distance < 1000 & Distance > -1000 & DonorCount > 1 & CTCFmotif != "noChange" & Class == "Other"][grep("MIR",Gene, invert = T)]
# mut2tss.counts[Distance < 1000 & Distance > -1000 & DonorCount > 1 & CTCFmotif != "noChange" & Class == "Other"][grep("MIR",Gene, invert = T)]


tmpdata <- mut2tss.counts[Distance < 1000 & Distance > -1000 & 
                              Gene %in% c("PAX5", "ZCCHC7", "MYC", "PVT1", "TERT","CLPTM1L", "NUCKS1","RAB29")]
tmpdata$Gene <- factor(tmpdata$Gene, levels = c("PAX5", "ZCCHC7", "MYC", "PVT1", "TERT","CLPTM1L", "NUCKS1","RAB29"))
tmpdata[is.na(CTCF1), CTCF1 := "noStrongMotifChange"]
tmpdata$CTCF1 <- factor(tmpdata$CTCF1, levels = levels(mut2tss.counts$CTCFmotif1))

tx <- tmpdata[Gene %in% c("CLPTM1L")]
tx[CTCF1 %in% c("noStrongMotifChange", "Gain"), CTCF2 := "noLoss"]
tx[CTCF1 %in% c("Loss"), CTCF2 := "Loss"]
tx$CTCF2 <- factor(tx$CTCF2, levels = c("noLoss", "Loss"))
p5 <- ggscatter(data = tx,
                x = "Distance", y = "DonorCount",
                color = "CTCF2", fill = "CTCF2", palette = get_palette(c("dodgerblue", "forestgreen"), 2),
                size = "ProjectCount",
                label = "Mutation", repel = T,
                font.label = list(color ="CTCF2"),
                label.select = list(criteria = "`CTCF2` == 'Loss'"),
                xlab = "Distance of Mutations to TSS", ylab = "Donor Counts",
                title = "-1000~1000, CLPTM1L, CTCF-motif-disrupting labeled") +
    scale_size(range = c(1, 2)) +
    geom_vline(xintercept = c(0), linetype = 1)
ggsave(filename = "ggscatter.ICGC.release28.WGS.refGene.tss.promoter.1k.count.dis.selectedgenes.CLPTM1L.pdf",
       useDingbats = F, dpi = 300, 
       width = 8, height = 3, ggpar(p5, ylim = c(0, max(tmpdata[Gene %in% c("CLPTM1L")]$DonorCount))))
