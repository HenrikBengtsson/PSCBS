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
setMethodS3("plotTracks", "CBS", function(x, scatter=TRUE, pch=20, col="gray", cex=1, grid=FALSE, Clim="auto", xScale=1e-6, Clab="auto", ..., byIndex=FALSE, add=FALSE) {
  # To please R CMD check
  fit <- x;


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'Clim':
  if (identical(Clim, "auto")) { 
    signalType <- getSignalType(fit);
    Clim <- switch(signalType,
      "log2ratio" = c(-3,3),
      "ratio"     = c(0,6),
      NULL
    );
  }

  if (identical(Clab, "auto")) { 
    signalType <- getSignalType(fit);
    Clab <- switch(signalType,
      "log2ratio" = "log2 CN ratio",
      "ratio"     = "CN ratio",
      NULL
    );
  }

  # Argument 'fit':
  if (nbrOfChromosomes(fit) > 1) {
    res <- plotTracksManyChromosomes(fit, scatter=scatter, pch=pch, Clim=Clim, xScale=xScale, Clab=Clab, ..., byIndex=byIndex, add=add);
    return(res);
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

  plot(x, CT, pch=pchT, cex=cex, col=col, ..., ylim=Clim, ylab=Clab);
  stext(side=3, pos=1, chrTag);
  if (grid) {
    yrange <- par("usr")[3:4];
    yrange[1] <- floor(yrange[1]);
    yrange[2] <- ceiling(yrange[2]);
    abline(h=seq(from=yrange[1], to=yrange[2], by=2), lty=3, col="gray");
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


setMethodS3("highlightCalls", "CBS", function(fit, pch=20, callCols=c(loss="red", gain="green", "amplification"="blue"), lwd=3, meanCol="purple", ..., xScale=1e-6, byIndex=FALSE, verbose=FALSE) {
  segs <- getSegments(fit, splitter=FALSE);

  # Identify segment calls
  callFields <- grep("Call$", colnames(segs));
  callTypes <- gsub("Call$", "", colnames(segs)[callFields]);
  nbrOfCalls <- length(callFields);

  # Nothing todo?
  if (nbrOfCalls == 0) {
    return();
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Textual and threshold annotations
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  params <- fit$params$callGainsAndLosses;
  if (!is.null(params)) {
    sigma <- params$sigmaMAD;
    scale <- params$scale;
    muR <- params$muR;

    # Draw thresholds used for calling
    abline(h=params$muR, col="gray", lty=3);
    abline(h=params$tauLoss, col=callCols["loss"], lty=3);
    abline(h=params$tauGain, col=callCols["gain"], lty=3);
    paramsGL <- substitute(paste(ng, " gains & ", nl, " losses [", hat(sigma)[MAD]==x, ", ", scale==y, ", ", hat(mu)[R]==z, "]"), list(ng=sum(segs$gainCall, na.rm=TRUE), nl=sum(segs$lossCall, na.rm=TRUE), x=round(sigma, digits=3), y=scale, z=muR));
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Tile chromosomes
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  fitT <- tileChromosomes(fit);
  verbose && str(verbose, fitT);



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Highlight gains and losses
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  dataT <- getLocusData(fitT);
  segsT <- getSegments(fitT, splitter=FALSE);
  chr <- dataT[,"chromosome"];
  x <- dataT[,"x"];
  y <- dataT[,3];
  nbrOfLoci <- nbrOfLoci(fitT);
  nbrOfSegments <- nbrOfSegments(fitT);
  rm(dataT);

  # For each non-neutral segment
  for (ss in seq(length=nbrOfSegments)) {
    seg <- segsT[ss,];

    for (tt in seq(along=callTypes)) {
      field <- callFields[tt];
      type <- callTypes[tt];

      # Called?
      call <- seg[[field]];
      if (isTRUE(call)) {
        col <- callCols[type];
        idxs <- which(chr == seg$chromosome & seg$start <= x & x <= seg$end);
        idxs <- Arguments$getIndices(idxs, max=nbrOfLoci);

        if (byIndex) {
          xs <- idxs;
        } else {
          xs <- x[idxs] * xScale;
        }
        ys <- y[idxs];
        points(xs, ys, pch=pch, col=col, ...);
        xx <- range(xs, na.rm=TRUE);
        yy <- rep(seg$mean, times=2);
        lines(xx, yy, lwd=lwd, col=meanCol);
      }
    } # for (tt ...)
  } # for (ss ...)
  rm(segsT);
}) # highlightCalls()


setMethodS3("getChromosomeOffsets", "CBS", function(fit, resolution=1e6, ...) {
  # Argument 'resolution':
  if (!is.null(resolution)) {
    resolution <- Arguments$getDouble(resolution, range=c(1,Inf));
  }

  data <- getChromosomeRanges(fit, ...);
  splits <- data[,"start"] + data[,"length"];

  if (!is.null(resolution)) {
    splits <- ceiling(splits / resolution);
    splits <- resolution * splits;
  }

  offsets <- c(0L, cumsum(splits));
  names(offsets) <- c(rownames(data), NA);

  offsets;
}, protected=TRUE) # getChromosomeOffsets()


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
  chrOffsets <- getChromosomeOffsets(fit, ...);
  chrStats <- getChromosomeRanges(fit);
  chrStats <- rbind(chrStats, c(1,NA,NA));
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

    verbose && exit(verbose);
  } # for (kk ...)

  # Update results
  fitT <- fit;
  fitT$data <- data;
  fitT$output <- segs;
  fitT$chromosomeStats <- chrStats;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Sanity checks
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  stopifnot(nbrOfLoci(fitT) == nbrOfLoci(fit));
  stopifnot(nbrOfSegments(fitT) == nbrOfSegments(fit));

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

  verbose && exit(verbose);

  fitT;
}, private=TRUE) # tileChromosomes()



setMethodS3("plotTracksManyChromosomes", "CBS", function(x, scatter=TRUE, pch=20, col="gray", Clim=c(0,6), xScale=1e-6, Clab="TCN", ..., subset=NULL, byIndex=FALSE, add=FALSE, onBegin=NULL, onEnd=NULL, verbose=FALSE) {
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
  fit <- tileChromosomes(fit);
  verbose && str(verbose, fit);
  # Sanity check
  stopifnot(!is.null(fit$chromosomeStats));

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
  chrStats <- fit$chromosomeStats;
  vs <- xScale * chrStats[,1:2];
  mids <- (vs[,1]+vs[,2])/2;
  CT <- y;

  nbrOfLoci <- length(x);
  chromosomes <- getChromosomes(fit);
  chrLabels <- sprintf("%02d", chromosomes);

  if (byIndex) {
    xs <- seq(along=x);
  } else {
    xs <- x;
  }

  if (!add) {
    par(mar=c(1,4,1,2)+1);
  }

  gh <- fit;
  gh$xScale <- xScale;

  xlim <- xScale*range(chrStats[,c("start","end")], na.rm=TRUE);
  xlab <- "Genomic position";

  pchT <- if (scatter) { pch } else { NA };

  plot(NA, xlim=xlim, ylim=Clim, xlab=xlab, ylab=Clab, axes=FALSE);
  if (!is.null(onBegin)) onBegin(gh=gh);
  points(xs, CT, pch=pchT, col=col, ...);
  mtext(text=chrLabels, side=rep(c(1,3), length.out=length(chrLabels)), at=mids, line=0.1, cex=0.9);
  abline(v=vs, lty=3);
  axis(side=2); box();
  drawLevels(fit, xScale=xScale, byIndex=byIndex);
  if (!is.null(onEnd)) onEnd(gh=gh);

  invisible(gh);
}, private=TRUE) # plotTracksManyChromosomes()


setMethodS3("drawCentromeres", "CBS", function(fit, genomeData, xScale=1e-6, col="gray", lty=3, ..., byIndex=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'genomeData':
  stopifnot(inherits(genomeData, "data.frame"));
  stopifnot(is.element("chrom", colnames(genomeData)));
  stopifnot(is.element("centroStart", colnames(genomeData)));
  stopifnot(is.element("centroEnd", colnames(genomeData)));

  # Calculate the midpoints of the centromeres
  genomeData$centroMid <- (genomeData$centroStart + genomeData$centroEnd) / 2;


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
  offsetsT <- offsets[seq(along=genomeData$centroMid)];
  midsT <- genomeData$centroMid + offsetsT;

  ats <- xScale*midsT;
  abline(v=ats, col=col, lty=lty, ...);

  invisible(ats);
}) # drawCentromeres()


############################################################################
# HISTORY:
# 2011-09-06
# o Added getChromosomeOffsets().
# 2011-09-01
# o BUG FIX: plotTracksManyChromosomes() for CBS gave an error because
#   internal variable 'CT' was not defined.
# o BUG FIX: tileChromosomes() for CBS not identify the chromosomes of
#   the loci, and hence generated corrupt/missing values while tiling.
# 2010-11-19
# o Created from PairedPSCBS.R.
############################################################################ 
