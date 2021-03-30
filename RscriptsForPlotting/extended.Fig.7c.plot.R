x <- read.table("otherlist/RAD21.5000.MCF-7.all.pairs.CTCF.FOXA1")

y <- x[x[,1]==x[,5],]
z <- data.frame(as.numeric(y[,2] < y[,6]),y[,c(11,14,17,20)])
colnames(z) <- c("EP_position","FOXA1.E","FOXA1.P","CTCF.E","CTCF.P")

library(reshape2)
library(ggplot2)

zzz <- z

zzz$FOXA1.E <- as.numeric(z$FOXA1.E)
zzz$FOXA1.P <- as.numeric(z$FOXA1.P)
zzz$CTCF.E <- as.numeric(z$CTCF.E)
zzz$CTCF.P <- as.numeric(z$CTCF.P)

zzz[,2:5] <- zzz[,2:5]-2


stats <- matrix(nrow=7,ncol=4)
rownames(stats) <- c("-><-","<-->","->->","<-<-","e.only","p.only","none")
colnames(stats) <- c("CTCF.E+CTCF.P", "FOXA1.E+CTCF.P", "FOXA1.E+FOXA1.P", "CTCF.E+FOXA1.P")

stats[1,1] <- nrow(zzz[zzz[,1]==1 & zzz[,4]==1 & zzz[,5]==-1,]) + nrow(zzz[zzz[,1]==0 & zzz[,4]==-1 & zzz[,5]==1,])
stats[1,2] <- nrow(zzz[zzz[,1]==1 & zzz[,2]==1 & zzz[,5]==-1,]) + nrow(zzz[zzz[,1]==0 & zzz[,2]==-1 & zzz[,5]==1,])
stats[1,3] <- nrow(zzz[zzz[,1]==1 & zzz[,2]==1 & zzz[,3]==-1,]) + nrow(zzz[zzz[,1]==0 & zzz[,2]==-1 & zzz[,3]==1,])
stats[1,4] <- nrow(zzz[zzz[,1]==1 & zzz[,4]==1 & zzz[,3]==-1,]) + nrow(zzz[zzz[,1]==0 & zzz[,4]==-1 & zzz[,3]==1,])

stats[2,1] <- nrow(zzz[zzz[,1]==1 & zzz[,4]==-1 & zzz[,5]==1,]) + nrow(zzz[zzz[,1]==0 & zzz[,4]==1 & zzz[,5]==-1,])
stats[2,2] <- nrow(zzz[zzz[,1]==1 & zzz[,2]==-1 & zzz[,5]==1,]) + nrow(zzz[zzz[,1]==0 & zzz[,2]==1 & zzz[,5]==-1,])
stats[2,3] <- nrow(zzz[zzz[,1]==1 & zzz[,2]==-1 & zzz[,3]==1,]) + nrow(zzz[zzz[,1]==0 & zzz[,2]==1 & zzz[,3]==-1,])
stats[2,4] <- nrow(zzz[zzz[,1]==1 & zzz[,4]==-1 & zzz[,3]==1,]) + nrow(zzz[zzz[,1]==0 & zzz[,4]==1 & zzz[,3]==-1,])

stats[3,1] <- nrow(zzz[zzz[,1]==1 & zzz[,4]==1 & zzz[,5]==1,]) + nrow(zzz[zzz[,1]==0 & zzz[,4]==-1 & zzz[,5]==-1,])
stats[3,2] <- nrow(zzz[zzz[,1]==1 & zzz[,2]==1 & zzz[,5]==1,]) + nrow(zzz[zzz[,1]==0 & zzz[,2]==-1 & zzz[,5]==-1,])
stats[3,3] <- nrow(zzz[zzz[,1]==1 & zzz[,2]==1 & zzz[,3]==1,]) + nrow(zzz[zzz[,1]==0 & zzz[,2]==-1 & zzz[,3]==-1,])
stats[3,4] <- nrow(zzz[zzz[,1]==1 & zzz[,4]==1 & zzz[,3]==1,]) + nrow(zzz[zzz[,1]==0 & zzz[,4]==-1 & zzz[,3]==-1,])

stats[4,1] <- nrow(zzz[zzz[,1]==1 & zzz[,4]==-1 & zzz[,5]==-1,]) + nrow(zzz[zzz[,1]==0 & zzz[,4]==1 & zzz[,5]==1,])
stats[4,2] <- nrow(zzz[zzz[,1]==1 & zzz[,2]==-1 & zzz[,5]==-1,]) + nrow(zzz[zzz[,1]==0 & zzz[,2]==1 & zzz[,5]==1,])
stats[4,3] <- nrow(zzz[zzz[,1]==1 & zzz[,2]==-1 & zzz[,3]==-1,]) + nrow(zzz[zzz[,1]==0 & zzz[,2]==1 & zzz[,3]==1,])
stats[4,4] <- nrow(zzz[zzz[,1]==1 & zzz[,4]==-1 & zzz[,3]==-1,]) + nrow(zzz[zzz[,1]==0 & zzz[,4]==1 & zzz[,3]==1,])

stats[5,1] <- nrow(zzz[zzz[,4]!=0 & zzz[,5]==0,])
stats[5,2] <- nrow(zzz[zzz[,2]!=0 & zzz[,5]==0,])
stats[5,3] <- nrow(zzz[zzz[,2]!=0 & zzz[,3]==0,])
stats[5,4] <- nrow(zzz[zzz[,4]!=0 & zzz[,3]==0,])

stats[6,1] <- nrow(zzz[zzz[,4]==0 & zzz[,5]!=0,])
stats[6,2] <- nrow(zzz[zzz[,2]==0 & zzz[,5]!=0,])
stats[6,3] <- nrow(zzz[zzz[,2]==0 & zzz[,3]!=0,])
stats[6,4] <- nrow(zzz[zzz[,4]==0 & zzz[,3]!=0,])

stats[7,1] <- nrow(zzz[zzz[,4]==0 & zzz[,5]==0,])
stats[7,2] <- nrow(zzz[zzz[,2]==0 & zzz[,5]==0,])
stats[7,3] <- nrow(zzz[zzz[,2]==0 & zzz[,3]==0,])
stats[7,4] <- nrow(zzz[zzz[,4]==0 & zzz[,3]==0,])

data <- melt(stats[1:4,])
pdf(file="ChIA-PET_paired_motif_pairs.orientation.pdf", width=6, height=6)
ggplot(data, aes(x = Var2, y = value, fill = factor(Var1))) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +
  xlab("") + ylab("pair count")+
  theme(legend.title = element_blank(), legend.position = "top")
dev.off()

data <- melt(stats)
pdf(file="ChIA-PET_paired_motif_pairs.orientation.all.pdf", width=6, height=6)
ggplot(data, aes(x = Var2, y = value, fill = factor(Var1))) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +
  xlab("") + ylab("pair count")+
  theme(legend.title = element_blank(), legend.position = "top")
dev.off()

