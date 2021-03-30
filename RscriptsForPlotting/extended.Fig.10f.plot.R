# screen -r tmp
# setwd("/data/shaojf/project.ERR/motif.mut")
library(data.table)
library(ggpubr)
# motif mutations ----------------------------------------------------------
cancersites <- fread("ICGC/ICGC.cancers", header = F)
setorder(cancersites, V3, V1)
motifmut2tss <- fread("ICGC.release28.WGS.refGene.tss.dis.sim.motif.txt", header = F)
motifmut2tss[, V4 := as.numeric(V4)]
motifmut2tss$V1 <- factor(motifmut2tss$V1, levels = cancersites$V1)

p1 <- gghistogram(data = motifmut2tss[abs(V4) < 5000, c(3,4,5)], 
                  x = "V4", y = "..count..", fill = "V5", color = "V5",
                  xlim = c(-5000,5000), bins = 100,
                  title = "all mutations (multiple specicimen will be counted repeatedly)")
ggsave(filename = "gghistogram.ICGC.release28.WGS.mut.motif.pdf", width = 15, height = 5, 
       ggpar(facet(p1, facet.by = "V5", nrow = 1), 
             xlab = "Distance of mutation to TSS", 
             ylab = "count of mutation"))


