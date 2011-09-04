###########################################################################/**
# @set "class=CBS"
# @RdocMethod plotTracks
#
# @title "Plots copy numbers along the genome"
#
# \description{
#  @get "title" for one or more chromosomes.
#  Each type of track is plotted in its own panel.
# }
#
# @synopsis
#
# \arguments{
#   \item{x}{A result object returned by @see "segmentByCBS".}
#   \item{pch}{The type of points to use.}
#   \item{Clim}{The range of copy numbers.}
#   \item{xScale}{The scale factor used for genomic positions.}
#   \item{...}{Not used.}
#   \item{add}{If @TRUE, the panels plotted are added to the existing plot,
#     otherwise a new plot is created.}
# }
#
# \value{
#   Returns nothing.
# }
# 
# @author
#
# @keyword IO
# @keyword internal
#*/########################################################################### 
setMethodS3("plotTracks", "CBS", function(x, scatter=TRUE, pch=20, col="gray", cex=1, grid=FALSE, Clim=c(0,6), xScale=1e-6, ..., add=FALSE) {
  # To please R CMD check
  fit <- x;
 
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'fit':
  if (nbrOfChromosomes(fit) > 1) {
    return(plotTracksManyChromosomes(fit, scatter=scatter, pch=pch, Clim=Clim, xScale=xScale, ..., add=add));
  }

  # Argument 'xScale':
  xScale <- Arguments$getNumeric(xScale, range=c(0,Inf));

  # Extract the input data
  data <- getLocusData(fit);
  if (is.null(data)) {
    throw("Cannot plot segmentation results. No input data available.");
  }

  chromosomes <- getChromosomes(fit);
  chromosome <- chromosomes[1];
  x <- data$x;
  CT <- data[,3];
  nbrOfLoci <- length(x);

  # Extract the segmentation
  segs <- fit$output;


  if (chromosome != 0) {
    chrTag <- sprintf("Chr%02d", chromosome);
  } else {
    chrTag <- "";
  }

  if (xScale != 1) {
    x <- xScale * x;
  }

  if (!add) {
    par(mar=c(1,4,1,2)+1);
  }

  pchT <- if (scatter) { pch } else { NA };

  plot(x, CT, pch=pchT, cex=cex, col=col, ylim=Clim, ylab="TCN");
  stext(side=3, pos=1, chrTag);
  if (grid) {
    abline(h=seq(from=0, to=Clim[2], by=2), lty=3, col="gray");
    abline(h=0, lty=1, col="black");
  }
  drawLevels(fit, col="purple", xScale=xScale);

  invisible();  
}) # plotTracks()


setMethodS3("plot", "CBS", function(x, ...) {
  plotTracks(x, ...);
})


setMethodS3("drawLevels", "CBS", function(fit, xScale=1e-6, ...) {
  # Get segmentation results
  segs <- as.data.frame(fit);

  # Extract subset of segments
  fields <- c("start", "end", "mean");
  segs <- segs[,fields, drop=FALSE];
  segs <- unique(segs);

  # Reuse drawLevels() for the DNAcopy class
  colnames(segs) <- c("loc.start", "loc.end", "seg.mean");
  dummy <- list(output=segs);
  class(dummy) <- "DNAcopy";
  drawLevels(dummy, xScale=xScale, ...);
})



setMethodS3("tileChromosomes", "CBS", function(fit, chrStarts=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'chrStarts':
  if (!is.null(chrStarts)) {
    chrStarts <- Arguments$getDoubles(chrStarts);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Tile chromosomes");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract data and segments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  data <- fit$data;
  segs <- fit$output;

  # Identify all chromosome
  chromosomes <- getChromosomes(fit);
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Additional chromosome annotations
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (is.null(chrStarts)) {
    xRange <- matrix(0, nrow=length(chromosomes), ncol=2);
    for (kk in seq(along=chromosomes)) {
      chromosome <- chromosomes[kk];
      idxs <- which(data$chromosome == chromosome);
      x <- data$x[idxs];
      r <- range(x, na.rm=TRUE);
      r <- r / 1e6;
      r[1] <- floor(r[1]);
      r[2] <- ceiling(r[2]);
      r <- 1e6 * r;
      xRange[kk,] <- r;
    } # for (kk ...)

    chrLength <- xRange[,2];
    chrStarts <- c(0, cumsum(chrLength)[-length(chrLength)]);
    chrEnds <- chrStarts + chrLength;

    # Not needed anymore
    rm(x, idxs);
  } # if (is.null(chrStarts))

  verbose && cat(verbose, "Chromosome starts:");
  chromosomeStats <- cbind(start=chrStarts, end=chrEnds, length=chrEnds-chrStarts);
  verbose && print(chromosomeStats);


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
    offset <- chrStarts[kk];
    verbose && cat(verbose, "Offset: ", offset);

    # Offset data
    idxs <- which(data$chromosome == chromosome);
    data$x[idxs] <- offset + data$x[idxs];

    # Offset segmentation
    idxs <- which(segs$chromosome == chromosome);
    segs[idxs,segFields] <- offset + segs[idxs,segFields];

    verbose && exit(verbose);
  } # for (kk ...)

  # Update results
  fit$data <- data;
  fit$output <- segs;
  fit$chromosomeStats <- chromosomeStats;

  verbose && exit(verbose);

  fit;
}) # tileChromosomes()



setMethodS3("plotTracksManyChromosomes", "CBS", function(x, pch=".", Clim=c(0,6), xScale=1e-6, ..., subset=NULL, add=FALSE, onBegin=NULL, onEnd=NULL, verbose=FALSE) {
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
  fit <- tileChromosomes(fit, verbose=verbose);
  verbose && str(verbose, fit);

  # Extract the input data
  data <- fit$data;
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
  vs <- xScale * fit$chromosomeStats[,1:2];
  mids <- (vs[,1]+vs[,2])/2;
  CT <- y;

  nbrOfLoci <- length(x);
  chromosomes <- getChromosomes(fit);
  chrTags <- sprintf("Chr%02d", chromosomes);

  if (!add) {
    par(mar=c(1,4,1,2)+1);
  }

  gh <- fit;
  gh$xScale <- xScale;

  xlim <- range(x, na.rm=TRUE);
  xlab <- "Genomic position";

  plot(NA, xlim=xlim, ylim=Clim, xlab=xlab, ylab="TCN", axes=FALSE);
  if (!is.null(onBegin)) onBegin(gh=gh);
  points(x, CT, pch=pch, col="gray");
  mtext(text=chrTags, side=rep(c(1,3), length.out=length(chrTags)), at=mids, line=0.1, cex=0.7);
  abline(v=vs, lty=3);
  axis(side=2); box();
  drawLevels(fit, xScale=xScale);
  if (!is.null(onEnd)) onEnd(gh=gh);

  invisible(gh);
}) # plotTracksManyChromosomes()



############################################################################
# HISTORY:
# 2011-09-01
# o BUG FIX: plotTracksManyChromosomes() for CBS gave an error because
#   internal variable 'CT' was not defined.
# o BUG FIX: tileChromosomes() for CBS not identify the chromosomes of
#   the loci, and hence generated corrupt/missing values while tiling.
# 2010-11-19
# o Created from PairedPSCBS.R.
############################################################################ 
