# screen -r tmp
# setwd("/data/shaojf/project.ERR/promoter.mut")
library(data.table)
library(ggpubr)
# read in general data ----------------------------------------------------------
cancersites <- fread("ICGC/ICGC.cancers", header = F)
setorder(cancersites, V3, V1)
mut2tss <- fread("ICGC/ICGC.release28.WGS.refGene.tss.dis.sim.txt", header = F)
mutinfo.1 <- fread("ICGC/ICGC.release28.WGS.mut.donor.count.1", header = T)
mutinfo.2 <- fread("ICGC/ICGC.release28.WGS.mut.donor.count.2", header = T)
mutinfo.1[, CancerDonors := mutinfo.2[match(mutinfo.1$Cancer, Cancer)]$DonorCount]
mutinfo.1[, SamplePerc := DonorCount/CancerDonors * 100]
mutcounts <- as.data.table(table(mutinfo.1$Cancer))
mutinfo.2[, UniqueMutationCount := mutcounts[match(mutinfo.2$Cancer,V1)]$N]
mutcounts1 <- as.data.table(table(unique(mut2tss[abs(V4) < 5000, c(1,3)])$V1))
mutinfo.2[, Unique5KPMutationCount := mutcounts1[match(mutinfo.2$Cancer,V1)]$N]
mutinfo.2[, Labs := paste0(Unique5KPMutationCount, "/", UniqueMutationCount)]

mut2tss[, Labs := mutinfo.2[match(mut2tss$V1, Cancer)]$Labs]
mut2tss[, Tissues := cancersites[match(mut2tss$V1, V1)]$V3]

plotdata <- mut2tss[V1 %in% mutinfo.2[DonorCount > 100 & UniqueMutationCount > 1000]$Cancer & abs(V4) < 5000]
plotdata.text <- mutinfo.2[DonorCount > 100 & UniqueMutationCount > 1000]
plotdata$V1 <- factor(plotdata$V1, levels = cancersites$V1)

p1 <- ggviolin(data = plotdata, 
               x = "V1", y = "V4", fill = "Tissues", color = "Tissues",
               ylim = c(-5000,5000), title = "mutations in cancers with DonorCount > 100 & UniqueMutationCount > 1000") + 
    geom_text(data = plotdata.text, aes(x = Cancer, y = 0, label = Labs), angle = 90)
ggsave(filename = "ggviolin.ICGC.release28.WGS.mut.stats.pdf", width = 12, height = 6, 
       ggpar(p1, x.text.angle = 30, xlab = "", ylab = "Distance of mutation to TSS"))

