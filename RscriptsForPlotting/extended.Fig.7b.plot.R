library(pheatmap)
x <- read.table("otherlist/RAD21.5000.MCF-7.all.pairs.final")

y <- x[x[,1]==x[,5],]
z <- data.frame(y[,c(4,8)],as.numeric(y[,2] < y[,6]),y[,c(11,14,17,20,23,26)])
colnames(z) <- c("E","P", "EP_position",
	"ESR1.E","ESR1.P","FOXA1.E","FOXA1.P","CTCF.E","CTCF.P")

zzz <- z
zzz$ESR1.E <- as.numeric(z$ESR1.E)
zzz$ESR1.P <- as.numeric(z$ESR1.P)
zzz$FOXA1.E <- as.numeric(z$FOXA1.E)
zzz$FOXA1.P <- as.numeric(z$FOXA1.P)
zzz$CTCF.E <- as.numeric(z$CTCF.E)
zzz$CTCF.P <- as.numeric(z$CTCF.P)

zzz[,4:9] <- zzz[,4:9]-2
a <- as.matrix(unique(zzz))[,4:9]
b <- matrix(as.numeric(a),ncol=6, byrow=F)
colnames(b) <- colnames(a)
b[b[,1]!=0,1] <- 1
b[b[,2]!=0,2] <- 2
b[b[,3]!=0,3] <- 3
b[b[,4]!=0,4] <- 4
b[b[,5]!=0,5] <- 5
b[b[,6]!=0,6] <- 6
pdf(file="ChIA-PET_paired_motif_pairs.pheatmap.pdf", width=5, height=5)
pheatmap(b[order(b[,6],b[,5],b[,3],b[,1],b[,4],b[,2],decreasing=T),c(1,3,5,6,4,2)], 
	color = colorRampPalette(c("grey","blue","red"))(7), cluster_rows = F, cluster_cols = F)
dev.off()
