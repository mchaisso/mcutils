library(getopt)
options <- matrix(c("dots", "d", 2, "character",
                    "xlabel", "x", 2, "character",
                    "ylabel", "y", 2, "character",
                    "colorize", "c", 0, "logical",
                    "width", "w", 2, "numeric",
                    "output", "o", 2, "character",
                    "querystart", "q", 2, "numeric",
                    "querytrack", "a", 2, "character",
                    "targetstart", "t", 2, "numeric",
                    "targettrack", "b", 2, "character"), byrow=T, ncol=4)

args <- getopt(options)
library(RColorBrewer)

pal=c("black", "red", "blue", "green")
t <- read.table(args$dots)

if (length(grep(".pdf",args$output)) > 0) {
  pdf(args$output)
} else if (length(grep(".png",args$output)) > 0) {
  options(bitmapType='cairo')
  png(args$output)
}
lineWidth = 0.01
if (is.null(args$width) == F) {
  lineWidth = args$width;
}
showAxes=T
if (is.null(args$querystart) == F) {
  showAxes=F
}


if (is.null(args$targetstart) == F) {
  showAxes=F
}

xRange <- range(t$V1)
yRange <- range(t$V2)
options(scipen=20)

plot(c(),xlim=xRange,ylim=yRange,xlab=args$xlabel,ylab=args$ylabel, axes=showAxes)
nCols <- dim(t)[2]

GetTicks <- function(start, span) {
  step = (10^max(2,floor(log10(span))))/2
  startTick <- floor(start/step)*step
  ticks <- c(start, seq(startTick + step, startTick+span, by=step))
  at <- c(ticks - start)
  labels <- paste(ticks)
  res <- data.frame(ticks, at,labels)
  return(res);
}
xSpan <- xRange[2] - xRange[1]
ySpan <- yRange[2] - yRange[1]

if (is.null(args$querystart) == F) {

  xTick <- GetTicks(args$querystart, xSpan)
  axis(1,at=xTick$at, labels=xTick$labels);
}
if (is.null(args$targetstart) == F) {
#  yTick <- GetTicks(args$targetstart, ySpan)
  yTick <- GetTicks(args$targetstart, ySpan)
#  yTick$at <- yTick$at+rep(min(t$V2),length(yTick$at))
  yTick$at <- yTick$at+min(t$V2)
  axis(2,at=yTick$at, labels=yTick$labels);
}

queryTrack <- c()
if (is.null(args$querytrack) == F) {
  if (file.info(args$querytrack)$size > 0) {
    queryTrack <- read.table(args$querytrack);
  }
}
targetTrack <- c()

if (is.null(args$targetrack) == F) {
  targetTrack <- read.table(args$targetTrack);
}


#
# Determine color by strand if supplied.
#
#useIndex = TRUE;

strands <- rep(0,length(t$V4))

for (idx in seq(min(t$V4),max(t$V4))) {
  i <- which(t$V4 == idx);
  lcol = pal[(idx%%4)+1 ];
  i = which(t$V4 == idx & t$V5 == 0) 	;
  segments(t$V1[i], t$V2[i], t$V1[i]+t$V3[i], t$V2[i]+t$V3[i], col=lcol, lwd=lineWidth);
  i <- which(t$V4 == idx & t$V5 == 1);
  segments(t$V1[i], t$V2[i], x1=t$V1[i]-t$V3[i], y1=t$V2[i]+t$V3[i], col=lcol, lwd=lineWidth);  
}

if (length(queryTrack) > 0) {
  par(xpd=NA)
  tHeight = ySpan*0.025
  nTrack <- length(queryTrack$V2)
  xB <- as.numeric(queryTrack$V2-args$querystart)
  xT <- as.numeric(queryTrack$V3-args$querystart)
  yB <- rep(tHeight*-0.25, nTrack)+min(t$V2)
  yT <- rep(tHeight*-1.25, nTrack)+min(t$V2)

#  lapply(seq(1,nTrack), function(i) rect(xB[i], yB[i], xT[i], yT[i], col="blue"))
  rect(xB,yB,xT,yT,col="blue")
  text(xB, yB, labels=queryTrack$V4, pos=2)
}


if (length(targetTrack) > 0) {

  par(xpd=NA)
  tHeight = ySpan*0.025

  nTrack <- length(queryTrack$V2)
  xB <- as.numeric(queryTrack$V2-args$querystart)
  xT <- as.numeric(queryTrack$V3-args$querystart)
  yB <- rep(tHeight*-0.25, nTrack)+min(t$V2)
  yT <- rep(tHeight*-1.25, nTrack)+min(t$V2)

#  lapply(seq(1,nTrack), function(i) rect(xB[i], yB[i], xT[i], yT[i], col="blue"))
  rect(xB,yB,xT,yT,col="blue")

}

dev.off()

outName <- paste(getwd(),"/",args$output,sep="")
print(outName)


