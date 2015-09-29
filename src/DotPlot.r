args <- commandArgs(trailingOnly = TRUE)

t <- read.table(args[1])
f <- which(t$V4 == 0 & t$V3 > 200)
r <- which(t$V4 == 1 & t$V3 > 200)

pdf(args[2])
plot(c(), xlim=c(min(t$V2),max(t$V2)), ylim=c(min(t$V1),max(t$V1)))
segments(t$V2[f], t$V1[f], t$V2[f]+t$V3[f], t$V1[f]+t$V3[f])
segments(t$V2[r], t$V1[r], t$V2[r]+t$V3[r], t$V1[r]+t$V3[r], col='red')
dev.off()
