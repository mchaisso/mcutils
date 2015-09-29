library(getopt)
options <- matrix(c("dots", "d", 2, "character",
                    "xlabel", "x", 2, "character",
                    "ylabel", "y", 2, "character",
                    "colorize", "c", 0, "logical",                    
                    "output", "o", 2, "character"), byrow=T, ncol=4)

args <- getopt(options)

pal=c("black", "red", "blue", "green")
t <- read.table(args$dots)

if (length(grep(".pdf",args$output)) > 0) {
  pdf(args$output)
} else if (length(grep(".png",args$output)) > 0) {
  png(args$output)
}

xRange <- range(t$V1)
yRange <- range(t$V2)
options(scipen=20)
plot(c(),xlim=xRange,ylim=yRange,xlab=args$xlabel,ylab=args$ylabel)
nCols <- dim(t)[2]
print(nCols)
if (nCols > 4) {
  strands <- t$V5;
  print(head(strands))
  print(tail(strands))
} else {
  strands <- rep(0,length(t$V4));
}
for (idx in seq(0,max(t$V4))) {
  i <- which(t$V4 == idx & strands == 0);
  print(head(strands == 0))
  color = pal[(idx%%4)+1 ];
  segments(t$V1[i], t$V2[i], t$V1[i]+t$V3[i], t$V2[i]+t$V3[i], col=color, lwd=0.01);
  i <- which(t$V4 == idx & strands == 1);
  print(head(i))
  segments(t$V1[i], t$V2[i], t$V1[i]-t$V3[i], t$V2[i]+t$V3[i], col=color, lwd=0.01);  
}

dev.off()

outName <- paste(getwd(),"/",args$output,sep="")
print(outName)


