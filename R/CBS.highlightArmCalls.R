setMethodS3("highlightArmCalls", "CBS", function(fit, genomeData, minFraction=0.95, callCols=c("loss"="red", "gain"="green"), ...) {
  # To please/trick R CMD check
  chromosome <- x <- NULL; rm(chromosome, x);

  # Argument 'minFraction':
  minFraction <- Arguments$getDouble(minFraction, range=c(0,1));

  regions <- getChromosomeRanges(fit);
  regions$end <- genomeData$centroStart;
  regions$start <- pmin(regions$start, regions$end);

  # Shrink regions
  for (rr in seq(length=nrow(regions))) {
    chr <- regions[rr,"chromosome"];
    x0 <- regions[rr,"start"];
    x1 <- regions[rr,"end"];
    xs <- subset(fit$data, chromosome == chr & x0 <= x & x <= x1)$x;
    if (length(xs) > 0) {
      range <- range(xs, na.rm=TRUE);
      x0 <- max(c(x0, range[1]), na.rm=TRUE);
      x1 <- min(c(x1, range[2]), na.rm=TRUE);
      regions[rr,"start"] <- x0;
      regions[rr,"end"] <- x1;
    }
  } # for (rr ...)
  regions[,"length"] <- regions[,"end"] - regions[,"start"] + 1L;
  callStatsP <- getCallStatistics(fit, regions=regions);
  

  regions <- getChromosomeRanges(fit);
  regions$start <- genomeData$centroEnd;
  regions$end <- pmax(regions$end, regions$start);

  # Shrink regions
  for (rr in seq(length=nrow(regions))) {
    chr <- regions[rr,"chromosome"];
    x0 <- regions[rr,"start"];
    x1 <- regions[rr,"end"];
    xs <- subset(fit$data, chromosome == chr & x0 <= x & x <= x1)$x;
    if (length(xs) > 0) {
      range <- range(xs, na.rm=TRUE);
      x0 <- max(c(x0, range[1]), na.rm=TRUE);
      x1 <- min(c(x1, range[2]), na.rm=TRUE);
      regions[rr,"start"] <- x0;
      regions[rr,"end"] <- x1;
    }
  } # for (rr ...)
  regions[,"length"] <- regions[,"end"] - regions[,"start"] + 1L;

  callStatsQ <- getCallStatistics(fit, regions=regions);

  callStats <- rbind(callStatsP, callStatsQ);

  # Not needed anymore
  rm(regions, regions, callStatsP, callStatsQ);

  keep <- c("chromosome", "start", "end", "lossFraction", "gainFraction");
  callStats <- callStats[,keep];
  colnames(callStats) <- gsub("Fraction", "", colnames(callStats));
  callTypes <- c("loss", "gain");

  # Adjust (start, end)
  offsets <- getChromosomeOffsets(fit);
  offsets <- offsets[callStats[,1]];
  callStats[,c("start","end")] <- offsets + callStats[,c("start","end")];

  usr <- par("usr");
  dy <- diff(usr[3:4]);
  xScale <- 1e-6;
  yy <- usr[3]+c(0,0.05*dy);
  abline(h=usr[3]+0.95*0.05*dy, lty=1, col="gray");

  for (key in callTypes) {
    xx <- callStats[,c("start", "end")];
    xx <- as.matrix(xx);
    xx <- xx * xScale;
    calls <- callStats[,key];
    for (kk in seq(along=calls)) {
      if (calls[kk] > 0) {
        lines(x=xx[kk,], y=rep(yy[1]+calls[kk]*0.05*dy, times=2), col=callCols[key]);
        if (calls[kk] > minFraction) {
          rect(xx[kk,1], yy[1], xx[kk,2], yy[2], col=callCols[key], border=NA);
        }
      }
    }
  } # for (key ...)
}, protected=TRUE); # highlightArmCalls()


############################################################################
# HISTORY:
# 2011-10-06
# o Created.
############################################################################
