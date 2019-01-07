library(getopt)
options <- matrix(c("dots", "d", 2, "character",
                    "xlabel", "x", 2, "character",
                    "ylabel", "y", 2, "character",
                    "colorize", "c", 0, "logical",
                    "width", "w", 2, "numeric",
                    "drawline", "l", 2, "logical",
                    "output", "o", 2, "character",
                    "querystart", "q", 2, "numeric",
                    "querytrack", "a", 2, "character",
                    "hline", "h", 2, "numeric",
                    "boxes", "", 2, "character",
                    "targetstart", "t", 2, "numeric",
                    "ystart", "e", 2, "numeric",
                    "yend", "z", 2, "numeric",
                    "xstart", "f", 2, "numeric",
                    "xend", "g", 2, "numeric",
                    "title", "1", 2, "character",
                    "targettrack", "b", 2, "character"), byrow=T, ncol=4)

args <- getopt(options)
library(RColorBrewer)

#pal=c("black", "red", "blue", "green")
nColors <- 9
pal <- c("black", brewer.pal(nColors,"Set1"))
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
  t$V1 = t$V1 + args$querystart
}
if (is.null(args$targetstart) == F) {
  t$V2 = t$V2 + args$targetstart
}

boxes <- c()
if (is.null(args$boxes) == F) {
  boxes <- read.table(args$boxes);
}

if (is.null(args$ystart) == F) {
  idx <- which(t$V2 >= args$ystart)
  t <- t[idx,]
}

if (is.null(args$yend) == F) {
  idx <- which(t$V2 <= args$yend)
  t <- t[idx,]
}

if (is.null(args$xstart) == F) {
  idx <- which(t$V1 >= args$xstart)
  t <- t[idx,]
}

if (is.null(args$xend) == F) {
  idx <- which(t$V1 <= args$xend)
  t <- t[idx,]
}


xRange <- range(t$V1)
xRange[1] <- max(0,xRange[1])
yRange <- range(t$V2)
yRange[1] <- max(0,yRange[1])
maxDim <- max(xRange[2],yRange[2])
xRange[2] <- maxDim
yRange[2] <- maxDim
options(scipen=20)
title <- ""
if (is.null(args$title) == F) {
  title <- args$title;
}
plot(c(),xlim=xRange,ylim=yRange,xlab=args$xlabel,ylab=args$ylabel, axes=showAxes, main=title)
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



queryTrack <- c()
if (is.null(args$querytrack) == F) {

  if (file.info(args$querytrack)$size > 0) {
    queryTrack <- read.table(args$querytrack, comment.char="");
  }
}
targetTrack <- c()

if (is.null(args$targettrack) == F) {
  if (file.info(args$targettrack)$size > 0) {  
    targetTrack <- read.table(args$targettrack, comment.char="");
  }
}


#
# Determine color by strand if supplied.
#
#useIndex = TRUE;

strands <- rep(0,length(t$V4))

for (idx in unique(t$V4)) {
  i <- which(t$V4 == idx);
  
  lcol = pal[(idx%%nColors)+1 ];
  i = which(t$V4 == idx & t$V5 == 0) 	;
  lineWidth <- 0.5

#  if (dim(t)[2] == 6) {
#    ltype=t$V6[i]
#  } else {
#    ltype=rep(1, length(i))
#  }
  ltype=1
  segments(t$V1[i], t$V2[i], t$V1[i]+t$V3[i], t$V2[i]+t$V3[i], col=lcol, lwd=lineWidth,lty=ltype);
  i <- which(t$V4 == idx & t$V5 == 1);
  segments(t$V1[i], t$V2[i], x1=t$V1[i]-t$V3[i], y1=t$V2[i]+t$V3[i], col=lcol, lwd=lineWidth,lty=ltype);

}

PlotGeneX <- function(b, y, height) {


  start <- as.integer(b$V2[1])
  end <- as.integer(b$V3[1])
  ##
  ## Plot gene models. Kind of a pain.
  ## 
  if (dim(b)[1] == 12) {
    sizes <- sapply(unlist(strsplit(as.character(b$V11[1]), ",")),as.integer)
    starts <- sapply(unlist(strsplit(as.character(b$V12[1]), ",")), as.integer) + start
  } else {
    starts <- c(start);
    sizes <- c(end-start);
  }
  if (dim(b)[1] == 1 & length(b) == 5){
    fillCol=as.character(b$V5[1]);
  } else {
    fillCol="grey";
  }
  rect( starts, y, starts+sizes, y+height, col=fillCol)
}


frac=0.025
if (length(queryTrack) > 0) {
  par(xpd=NA)
  tHeight = ySpan*frac

  nTrack <- length(queryTrack$V2)


  nRows <- 25
  n <- length(queryTrack$V2)
  yOffsets <- rep(0,n)
  lastX <- rep(0,n)
  for (j in seq(1,n)) {
    curRow <- 1
    while(curRow < n & lastX[curRow] > queryTrack$V2[j]) {
      curRow <- curRow+1;
    }
    lastX[curRow] <- queryTrack$V3[j];
    yOffsets[j] <- curRow*tHeight;
  }

  
  xB <- as.numeric(queryTrack$V2)
  xT <- as.numeric(queryTrack$V3)
  yB <- rep(tHeight*-0.25, nTrack)+min(t$V2) + yOffsets
  yT <- rep(tHeight*-1.25, nTrack)+min(t$V2) + yOffsets
  yBText <- rep(tHeight, nTrack)+min(t$V2)  
  h <- abs(yB[1]-yT[1])
  rect(xB,yT+h/2,xT,yT+h/2,col=t$V5)

  lapply(seq(1,length(xB)), function(i) PlotGeneX(queryTrack[i,], yT[i], h))
  lapply(seq(1,length(xB)), function(i) print(queryTrack[i,]))

  text(xB,yT +(yB-yT)*0.35, labels=queryTrack$V4, pos=2)


  
#  if (args$drawline == T) {
#    lapply(seq(1,length(yB)), function(i) {segments(xB[i],0,xB[i],yRange[2],lwd=0.1);})
#  }
  
}



PlotGeneY <- function(b, x, width, yOffset) {

  start <- as.integer(b$V2[1])
  end <- as.integer(b$V3[1])
  if (dim(b)[1]==12) {
    sizes <- sapply(unlist(strsplit(as.character(b$V11[1]), ",")),as.integer)
    starts <- sapply(unlist(strsplit(as.character(b$V12[1]), ",")), as.integer)
  } else {
    starts <- c(start)
    sizes <- c(end-start)
  }
  if (dim(b)[1] == 1 & length(b) == 5) {
    fillCol=as.character(b$V5[1]);
  } else {
    fillCol="gray";
  }
  rect(x, starts, x+width, starts+sizes, col=fillCol)
  text(x+2*width, end+yOffset, labels=b$V4[1], pos=2, srt=-90)  
}

if (length(targetTrack) > 0) {

#  par(xpd=NA)

  qWidth <- xSpan*frac

  ##
  ## Determine layering so that genes do not overlap.
  ##
  n <- length(targetTrack$V1)
  trackLastY <- rep(0, n)
  
  xOffset <- rep(0,n)
  xRange <- max(t$V2)-min(t$V2)
  yRange <- max(t$V3)-min(t$V3)
  yOffset <- yRange*0.01
  # Give some space for text
  textRange <- xRange*0.1
  for (j in seq(1,n)) {
    curTrack <- 1
    while (trackLastY[curTrack] > targetTrack$V2[j] & curTrack < n) {
      curTrack <- curTrack+1;
    }
    trackLastY[curTrack] <- targetTrack$V3[j]
    xOffset[j] <- qWidth*curTrack;
  }

  nTrack <- length(targetTrack$V2)
#  xOffsets <- xOffset[(seq(1,nTrack)%%(nCols+1))+1]

  xB <- rep(qWidth*0.25, nTrack)+min(t$V2) + xOffset
  xT <- rep(qWidth*1.25, nTrack)+min(t$V2) + xOffset
  yB <- as.numeric(targetTrack$V2)
  yT <- as.numeric(targetTrack$V3)
  xBText <- rep(qWidth, nTrack)+min(t$V1)
  w <- abs((xT-xB)/2)
  lapply(seq(1,length(targetTrack$V1)), function(tt) PlotGeneY(targetTrack[tt,], xT[tt], w, yOffset))
  

  if (is.null(args$drawline) == F) {
    idx = which(targetTrack$V4 !="")
    lapply(idx, function(i) { segments(xRange[1], yB[i], x1=xRange[2], lty=2,lwd=0.1)})
  }
}
if (is.null(args$hline) == F) {
  abline(h=args$hline,lwd=0.1);
}

lapply(seq(1,length(boxes)), function(i) rect(boxes$V1[i],boxes$V2[i], boxes$V3[i], boxes$V4[i], lwd=0.1))


dev.off()

outName <- paste(getwd(),"/",args$output,sep="")



