setMethodS3("tileChromosomes", "CBS", function(fit, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    pushState(verbose)
    on.exit(popState(verbose))
  }


  # Nothing to do, i.e. already tiled?
  if (isTRUE(attr(fit, "tiledChromosomes"))) {
    return(fit)
  }

  verbose && enter(verbose, "Tiling chromosomes")

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract data and segments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  data <- getLocusData(fit)
  segs <- getSegments(fit)
  knownSegments <- fit$params$knownSegments

  # Identify all chromosome
  chromosomes <- getChromosomes(fit)

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Additional chromosome annotations
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  chrStats <- getChromosomeRanges(fit)

  # Build an "empty" row with start == 1.
  chrStatsKK <- chrStats[1,,drop=FALSE][NA,,drop=FALSE]
  chrStatsKK[,"start"] <- 1L

  # Append empty row
  chrStats <- rbind(chrStats, chrStatsKK)

  # Offset (start,stop)
  chrOffsets <- getChromosomeOffsets(fit, ...)
  chrStats[,"start"] <- chrStats[,"start"] + chrOffsets
  chrStats[,"end"] <- chrStats[,"end"] + chrOffsets


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Offset...
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  segFields <- grep("(start|end)$", colnames(segs), value=TRUE)
  for (kk in seq_along(chromosomes)) {
    chromosome <- chromosomes[kk]
    chrTag <- sprintf("Chr%02d", chromosome)
    verbose && enter(verbose, sprintf("Chromosome #%d ('%s') of %d",
                                         kk, chrTag, length(chromosomes)))

    # Get offset for this chromosome
    offset <- chrOffsets[kk]
    verbose && cat(verbose, "Offset: ", offset)

    # Offset data
    idxs <- which(data$chromosome == chromosome)
    if (length(idxs) > 0L) {
      data$x[idxs] <- offset + data$x[idxs]
    }

    # Offset segmentation
    idxs <- which(segs$chromosome == chromosome)
    if (length(idxs) > 0L) {
      segs[idxs,segFields] <- offset + segs[idxs,segFields]
    }

    # Offset known segments
    idxs <- which(knownSegments$chromosome == chromosome)
    if (length(idxs) > 0L) {
      knownSegments[idxs,c("start", "end")] <- offset + knownSegments[idxs,c("start", "end")]
    }

    verbose && exit(verbose)
  } # for (kk ...)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Update results
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  fitT <- fit
  fitT$data <- data
  fitT$output <- segs
  fitT$chromosomeStats <- chrStats
  fitT$params$knownSegments <- knownSegments
  fitT$params$chrOffsets <- chrOffsets


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Sanity checks
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  segs <- getSegments(fit)
  segsT <- getSegments(fitT)
  segs <- segs[,!is.element(colnames(segs), c("start", "end"))]
  segsT <- segsT[,!is.element(colnames(segsT), c("start", "end"))]
  .stop_if_not(all.equal(segsT, segs))

  data <- getLocusData(fit)
  dataT <- getLocusData(fitT)
  data <- data[,!is.element(colnames(data), c("x"))]
  dataT <- dataT[,!is.element(colnames(dataT), c("x"))]
  .stop_if_not(all.equal(dataT, data))

  .stop_if_not(nbrOfLoci(fitT) == nbrOfLoci(fit))
  .stop_if_not(nbrOfSegments(fitT) == nbrOfSegments(fit))

  # Flag object
  attr(fitT, "tiledChromosomes") <- TRUE

  verbose && exit(verbose)

  fitT
}, protected=TRUE) # tileChromosomes()



setMethodS3("plotTracksManyChromosomes", "CBS", function(x, scatter=TRUE, pch=20, col="gray", meanCol="purple", Clim=c(0,3*ploidy(x)), xScale=1e-6, xlab="Genomic position", Clab="TCN", ..., boundaries=TRUE, levels=TRUE, subset=NULL, byIndex=FALSE, add=FALSE, onBegin=NULL, onEnd=NULL, mar=NULL, verbose=FALSE) {
  # To please R CMD check
  fit <- x

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'fit':

  # Argument 'add':
  add <- Arguments$getLogical(add)

  # Argument 'Clim':
  if (!add) {
    Clim <- Arguments$getNumerics(Clim, length=c(2L,2L),
                                        disallow=c("Inf", "NA", "NaN"))
  }

  # Argument 'xScale':
  xScale <- Arguments$getNumeric(xScale, range=c(0,Inf))

  # Argument 'subset':
  if (!is.null(subset)) {
    subset <- Arguments$getDouble(subset, range=c(0,1))
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    pushState(verbose)
    on.exit(popState(verbose))
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Tile chromosomes
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  fitT <- tileChromosomes(fit)
  verbose && str(verbose, fitT)
  # Sanity check
  .stop_if_not(!is.null(fitT$chromosomeStats))

  # Extract the input data
  data <- getLocusData(fitT)
  if (is.null(data)) {
    stop("Cannot plot segmentation results. No input data available.")
  }

  # Subset of the loci?
  if (!is.null(subset) && subset < 1) {
    n <- nrow(data)
    keep <- sample(n, size=subset*n)
    data <- data[keep,]
  }

  # To please R CMD check
  CT <- y <- muN <- betaT <- betaN <- betaTN <- NULL
  rm(list=c("CT", "muN", "betaT", "betaN", "betaTN"))
  attachLocally(data)
  x <- xScale * x
  chrStats <- fitT$chromosomeStats
  chrStats <- chrStats[-nrow(chrStats),,drop=FALSE]
  chrRanges <- as.matrix(chrStats[,c("start","end")])
  vs <- xScale * chrRanges
  mids <- (vs[,1]+vs[,2])/2
  CT <- y

  nbrOfLoci <- length(x)
  chromosomes <- getChromosomes(fitT)
  chrLabels <- sprintf("%02d", chromosomes)

  if (byIndex) {
    xs <- seq_along(x)
  } else {
    xs <- x
  }

  if (!add && !is.null(mar)) {
    par(mar=mar)
  }

  gh <- fitT
  gh$xScale <- xScale

  xlim <- xScale*range(chrRanges, na.rm=TRUE)

  pchT <- if (scatter) { pch } else { NA }

  plot(NA, xlim=xlim, ylim=Clim, xlab=xlab, ylab=Clab, axes=FALSE)
  if (!is.null(onBegin)) onBegin(gh=gh)
  points(xs, CT, pch=pchT, col=col, ...)
  side <- rep(c(1,3), length.out=length(chrLabels))
  mtext(text=chrLabels, side=side, at=mids, line=0.1, cex=0.7*par("cex"))
  if (boundaries) {
    abline(v=vs, lty=3)
  }
  axis(side=2)
  box()
  if (levels) {
    drawLevels(fitT, col=meanCol, xScale=xScale, byIndex=byIndex)
  }
  if (!is.null(onEnd)) onEnd(gh=gh)

  invisible(gh)
}, private=TRUE) # plotTracksManyChromosomes()
