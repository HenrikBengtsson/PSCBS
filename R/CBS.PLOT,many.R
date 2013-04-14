setMethodS3("tileChromosomes", "CBS", function(fit, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Tiling chromosomes");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract data and segments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  data <- getLocusData(fit);
  segs <- getSegments(fit);
  knownSegments <- fit$params$knownSegments;

  # Identify all chromosome
  chromosomes <- getChromosomes(fit);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Additional chromosome annotations
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  chrStats <- getChromosomeRanges(fit);

  # Build an "empty" row with start == 1.
  chrStatsKK <- chrStats[1,,drop=FALSE][NA,,drop=FALSE];
  chrStatsKK[,"start"] <- 1L;

  # Append empty row
  chrStats <- rbind(chrStats, chrStatsKK);

  # Offset (start,stop)
  chrOffsets <- getChromosomeOffsets(fit, ...);
  chrStats[,"start"] <- chrStats[,"start"] + chrOffsets;
  chrStats[,"end"] <- chrStats[,"end"] + chrOffsets;


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Offset...
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  segFields <- grep("(start|end)$", colnames(segs), value=TRUE);
  for (kk in seq(along=chromosomes)) {
    chromosome <- chromosomes[kk];
    chrTag <- sprintf("Chr%02d", chromosome);
    verbose && enter(verbose, sprintf("Chromosome #%d ('%s') of %d",
                                         kk, chrTag, length(chromosomes)));

    # Get offset for this chromosome
    offset <- chrOffsets[kk];
    verbose && cat(verbose, "Offset: ", offset);

    # Offset data
    idxs <- which(data$chromosome == chromosome);
    data$x[idxs] <- offset + data$x[idxs];

    # Offset segmentation
    idxs <- which(segs$chromosome == chromosome);
    segs[idxs,segFields] <- offset + segs[idxs,segFields];

    # Offset known segments
    idxs <- which(knownSegments$chromosome == chromosome);
    knownSegments[idxs,c("start", "end")] <- offset + knownSegments[idxs,c("start", "end")];

    verbose && exit(verbose);
  } # for (kk ...)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Update results
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  fitT <- fit;
  fitT$data <- data;
  fitT$output <- segs;
  fitT$chromosomeStats <- chrStats;
  fitT$params$knownSegments <- knownSegments;
  fitT$params$chrOffsets <- chrOffsets;


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Sanity checks
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  segs <- getSegments(fit);
  segsT <- getSegments(fitT);
  segs <- segs[,!is.element(colnames(segs), c("start", "end"))];
  segsT <- segsT[,!is.element(colnames(segsT), c("start", "end"))];
  stopifnot(all.equal(segsT, segs));

  data <- getLocusData(fit);
  dataT <- getLocusData(fitT);
  data <- segs[,!is.element(colnames(data), c("x"))];
  dataT <- segsT[,!is.element(colnames(dataT), c("x"))];
  stopifnot(all.equal(dataT, data));

  stopifnot(nbrOfLoci(fitT) == nbrOfLoci(fit));
  stopifnot(nbrOfSegments(fitT) == nbrOfSegments(fit));

  verbose && exit(verbose);

  fitT;
}, protected=TRUE) # tileChromosomes()



setMethodS3("plotTracksManyChromosomes", "CBS", function(x, scatter=TRUE, pch=20, col="gray", meanCol="purple", Clim=c(0,6), xScale=1e-6, xlab="Genomic position", Clab="TCN", ..., boundaries=TRUE, levels=TRUE, subset=NULL, byIndex=FALSE, add=FALSE, onBegin=NULL, onEnd=NULL, mar=NULL, verbose=FALSE) {
  # To please R CMD check
  fit <- x;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'fit':

  # Argument 'xScale':
  xScale <- Arguments$getNumeric(xScale, range=c(0,Inf));

  # Argument 'subset':
  if (!is.null(subset)) {
    subset <- Arguments$getDouble(subset, range=c(0,1));
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Tile chromosomes
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  fitT <- tileChromosomes(fit);
  verbose && str(verbose, fitT);
  # Sanity check
  stopifnot(!is.null(fitT$chromosomeStats));

  # Extract the input data
  data <- getLocusData(fitT);
  if (is.null(data)) {
    throw("Cannot plot segmentation results. No input data available.");
  }

  # Subset of the loci?
  if (!is.null(subset) && subset < 1) {
    n <- nrow(data);
    keep <- sample(n, size=subset*n);
    data <- data[keep,];
  }

  # To please R CMD check
  CT <- y <- muN <- betaT <- betaN <- betaTN <- NULL;
  rm(CT, muN, betaT, betaN, betaTN);
  attachLocally(data);
  x <- xScale * x;
  chrStats <- fitT$chromosomeStats;
  chrStats <- chrStats[-nrow(chrStats),,drop=FALSE];
  chrRanges <- as.matrix(chrStats[,c("start","end")]);
  vs <- xScale * chrRanges;
  mids <- (vs[,1]+vs[,2])/2;
  CT <- y;

  nbrOfLoci <- length(x);
  chromosomes <- getChromosomes(fitT);
  chrLabels <- sprintf("%02d", chromosomes);

  if (byIndex) {
    xs <- seq(along=x);
  } else {
    xs <- x;
  }

  if (!add && !is.null(mar)) {
    par(mar=mar);
  }

  gh <- fitT;
  gh$xScale <- xScale;

  xlim <- xScale*range(chrRanges, na.rm=TRUE);

  pchT <- if (scatter) { pch } else { NA };

  plot(NA, xlim=xlim, ylim=Clim, xlab=xlab, ylab=Clab, axes=FALSE);
  if (!is.null(onBegin)) onBegin(gh=gh);
  points(xs, CT, pch=pchT, col=col, ...);
  side <- rep(c(1,3), length.out=length(chrLabels));
  mtext(text=chrLabels, side=side, at=mids, line=0.1, cex=0.7*par("cex"));
  if (boundaries) {
    abline(v=vs, lty=3);
  }
  axis(side=2);
  box();
  if (levels) {
    drawLevels(fitT, col=meanCol, xScale=xScale, byIndex=byIndex);
  }
  if (!is.null(onEnd)) onEnd(gh=gh);

  invisible(gh);
}, private=TRUE) # plotTracksManyChromosomes()



setMethodS3("drawChromosomes", "CBS", function(x, lty=3, xScale=1e-6, ..., byIndex=FALSE, verbose=FALSE) {
  # To please R CMD check
  fit <- x;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'fit':

  # Argument 'xScale':
  xScale <- Arguments$getNumeric(xScale, range=c(0,Inf));

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Tile chromosomes
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  fitT <- tileChromosomes(fit);
  verbose && str(verbose, fitT);
  # Sanity check
  stopifnot(!is.null(fitT$chromosomeStats));

  chrStats <- fitT$chromosomeStats;
  chrStats <- chrStats[-nrow(chrStats),,drop=FALSE];
  chrRanges <- as.matrix(chrStats[,c("start","end")]);
  vs <- xScale * chrRanges;
  mids <- (vs[,1]+vs[,2])/2;
  chromosomes <- getChromosomes(fitT);
  chrLabels <- sprintf("%02d", chromosomes);
  side <- rep(c(1,3), length.out=length(chrLabels));
  mtext(text=chrLabels, side=side, at=mids, line=0.1, cex=0.7*par("cex"));
  abline(v=vs, lty=lty);
}, protected=TRUE) # drawChromosomes()



setMethodS3("drawCentromeres", "CBS", function(fit, genomeData, what=c("start", "end"), xScale=1e-6, col="gray", lty=3, ..., byIndex=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'genomeData':
  stopifnot(inherits(genomeData, "data.frame"));
  stopifnot(is.element("chromosome", colnames(genomeData)));
  stopifnot(is.element("centroStart", colnames(genomeData)));
  stopifnot(is.element("centroEnd", colnames(genomeData)));

  # Calculate the midpoints of the centromeres
  colnames(genomeData) <- tolower(gsub("centro", "", colnames(genomeData)));
  genomeData$mid <- (genomeData$start + genomeData$end) / 2;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Tile chromosomes
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  fitT <- tileChromosomes(fit);
  verbose && str(verbose, fitT);
  # Sanity check
  stopifnot(!is.null(fitT$chromosomeStats));

  chrStats <- fitT$chromosomeStats;
  offsets <- chrStats[,"start"] - chrStats[1,"start"];

  # Centroid locations in the tiled space
  offsetsT <- offsets[seq(length=nrow(genomeData))];

  xx <- genomeData[,what,drop=FALSE];
  xx <- as.matrix(xx);
  xx <- offsetsT + xx;

  ats <- xScale * xx;
  for (cc in seq(length=ncol(xx))) {
    abline(v=ats[,cc], col=col, lty=lty, ...);
  }

  invisible(ats);
}, protected=TRUE) # drawCentromeres()


setMethodS3("highlightArmCalls", "CBS", function(fit, genomeData, minFraction=0.95, callCols=c("loss"="red", "gain"="green"),   xScale=1e-6, ...) {
  # To please/trick R CMD check
  chromosome <- x <- NULL; rm(chromosome, x);

  callStats <- callArms(fit, genomeData=genomeData, minFraction=minFraction);

  callTypes <- grep("Fraction", colnames(callStats), value=TRUE);
  callTypes <- gsub("Fraction", "", callTypes);

  callTypes <- intersect(callTypes, c("loss", "gain"));

  # Adjust (start, end)
  offsets <- getChromosomeOffsets(fit);
  offsets <- offsets[callStats[,"chromosome"]];
  callStats[,c("start","end")] <- offsets + callStats[,c("start","end")];

  nbrOfRegions <- nrow(callStats);

  # Nothing todo?
  if (nbrOfRegions == 0) {
    return(invisible(callStats));
  }


  usr <- par("usr");
  dy <- diff(usr[3:4]);
  yy <- usr[3]+c(0,0.05*dy);
  abline(h=usr[3]+0.95*0.05*dy, lty=1, col="gray");

  xx <- callStats[,c("start", "end")];
  xx <- as.matrix(xx);
  xx <- xx * xScale;

  for (type in callTypes) {
    col <- callCols[type];
    keyA <- sprintf("%sFraction", type);
    keyB <- sprintf("%sCall", type);
    for (kk in seq(length=nbrOfRegions)) {
      xs <- xx[kk,];
      score <- callStats[kk, keyA];
      if (is.finite(score) && score > 0) {
        ys <- rep(yy[1]+callStats[kk, keyA]*0.05*dy, times=2);
        lines(x=xs, y=ys, col=col);
        call <- callStats[kk, keyB];
        if (call) {
          rect(xs[1], yy[1], xs[2], yy[2], col=col, border=NA);
        }
      }
    }
  } # for (type ...)

  invisible(callStats);
}, protected=TRUE); # highlightArmCalls()




############################################################################
# HISTORY:
# 2011-12-06
# o Now plotTracks() for CBS always returns an invisible object.
# 2011-12-03
# o Added drawChangePoints() for CBS.
# 2011-10-23
# o BUG FIX: highlightArmCalls() for CBS did not handle empty chromosomes.
# 2011-10-08
# o Added drawChromosomes() for CBS.
# 2011-10-07
# o Added highlightArmCalls() for CBS.
# 2011-10-06
# o Now drawCentromeres() for CBS can also plot start and stop.
# o ROBUSTNESS: Now plotTracksManyChromosomes() extract (start,end)
#   information by names (no longer assuming certain indices).
# 2011-09-07
# o Added highlightLocusCalls().
# 2011-09-06
# o Added highlightCalls().
# o Added getChromosomeOffsets().
# 2011-09-01
# o BUG FIX: plotTracksManyChromosomes() for CBS gave an error because
#   internal variable 'CT' was not defined.
# o BUG FIX: tileChromosomes() for CBS not identify the chromosomes of
#   the loci, and hence generated corrupt/missing values while tiling.
# 2010-11-19
# o Created from PairedPSCBS.R.
############################################################################
