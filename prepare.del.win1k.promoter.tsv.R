# setwd("/data/shaojf/project.ERR/promoter.del")
library(data.table)
del.win1k <- fread("ICGC.release28.del.refGene.promoter.txt", header = F, sep = "\t")
del.win1k[, DelLength := V3-V2]
del.win1k[, Class := ifelse(DelLength < 1e4, "pSmallDEL","pLargeDEL")]
del.win1k[, Project := matrix(unlist(strsplit(V4, split = "\\|")), byrow = T, ncol = 4)[,1]]
del.win1k[, Donor := matrix(unlist(strsplit(V4, split = "\\|")), byrow = T, ncol = 4)[,2]]
del.win1k[, Specimen := matrix(unlist(strsplit(V4, split = "\\|")), byrow = T, ncol = 4)[,3]]
del.win1k[, Gene := tss2gene[match(del.win1k$V8, V2)]$V1]
del.win1k.p <- del.win1k[grep("MIR",V8, invert = T)]
del.win1k.p[Gene %in% cosmics$`Gene Symbol`, CosmicGene := Gene]
del.win1k.p[, ClosetoCosmicGene := onconeighbors[match(del.win1k.p$Gene, V2)]$V1]
del.win1k.p[, CosmicGeneRole := cosmics[match(del.win1k.p$CosmicGene, `Gene Symbol`)]$`Role in Cancer`]
del.win1k.p[, ClosetoCosmicGeneRole := cosmics[match(del.win1k.p$ClosetoCosmicGene, `Gene Symbol`)]$`Role in Cancer`]
fwrite(x = del.win1k.p, file = "del.win1k.promoter.tsv", sep = "\t")

