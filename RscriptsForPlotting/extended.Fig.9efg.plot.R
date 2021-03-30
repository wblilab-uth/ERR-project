# setwd("/data/shaojf/project.ERR/promoter.mut")
library(maftools)
library(data.table)

#####
# final
# TERT	SLC12A7
# TERT	CTD-3080P12.3
# TERT	SLC6A19
# TERT	SLC6A18
# TERT	MIR4457
# TERT	CLPTM1L
# TERT	LINC01511
# TERT	SLC6A3
# MYC	CASC11
# MYC	PVT1
# MYC	MIR1204
# PAX5	MIR4540
# PAX5	MIR4476
# PAX5	EBLN3P
# PAX5	ZCCHC7
# PAX5	ZCCHC7

col <- c("firebrick4")
names(col) = c("TSS+(-1000,1000)")
mymaf <- read.maf(maf = "ICGC/ICGC.release28.WGS.refGene.tss.promoter.maf", vc_nonSyn = c("TSS+(-1000,1000)"))

pdf(file = paste("ICGC.release28.WGS.TERT", "oncoplot.1col.pdf", sep = "."), width = 10, height = 4)
oncoplot(maf = mymaf, colors = col, fontSize = 0.7,
         genes = c("TERT","CLPTM1L", "SLC12A7", "SLC6A19", "SLC6A18", "CTD-3080P12.3", "LINC01511", "SLC6A3"))
dev.off()

pdf(file = paste("ICGC.release28.WGS.MYC", "oncoplot.1col.pdf", sep = "."), width = 10, height = 2)
oncoplot(maf = mymaf, colors = col, fontSize = 0.7, genes = c("MYC", "PVT1", "CASC11"))
dev.off()

pdf(file = paste("ICGC.release28.WGS.PAX5", "oncoplot.1col.pdf", sep = "."), width = 10, height = 2)
oncoplot(maf = mymaf, colors = col, fontSize = 0.7, genes = c("PAX5", "ZCCHC7", "EBLN3P"))
dev.off()

