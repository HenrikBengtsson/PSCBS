#   \item{chromosomes}{An optional @numeric @vector specifying which
#     chromosomes to plot.}
#
#   \item{seed}{An (optional) @integer specifying the random seed to be
#     set before subsampling.  The random seed is
#     set to its original state when exiting.  If @NULL, it is not set.}
#
#   \item{verbose}{See @see "R.utils::Verbose".}
#
setMethodS3("plotTracksManyChromosomes", "PairedPSCBS", function(fit, chromosomes=getChromosomes(fit), tracks=NULL, scatter="*", calls=if (callLoci || length(chromosomes) == 1L) ".*" else NULL, callLoci=FALSE, callThresholds=TRUE, boundaries=TRUE, knownSegments=FALSE, quantiles=c(0.05,0.95), seed=0xBEEF, pch=".", Clim=c(0,3*ploidy(fit)), Blim=c(0,1), xScale=1e-6, xlabTicks=if (length(chromosomes) == 1L) "[pos]" else "[chr]", ..., subset=if (length(chromosomes) > 1L) 0.1 else NULL, add=FALSE, subplots=!add && (length(tracks) > 1L), oma=c(0,0,2,0), mar=c(2,5,1,3)+0.1, onBegin=NULL, onEnd=NULL, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  attachGH <- function(gh, envir=parent.frame()) {
    if (!is.list(gh)) return()
    if (is.null(gh$track)) return()
    if (!is.null(value <- gh$track)) assign("track", value, envir=envir)
    if (!is.null(value <- gh$subtracks)) assign("trackT", value, envir=envir)
    if (!is.null(value <- gh$scatter$col)) assign("colS", value, envir=envir)
    if (!is.null(value <- gh$scatter$pch)) assign("pchT", value, envir=envir)
    if (!is.null(value <- gh$level$col)) assign("colL", value, envir=envir)
    if (!is.null(value <- gh$cis$col)) assign("colC", value, envir=envir)
  } # attachGH()


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Graphical styles
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  opts <- list(
    scatter = list(pch=".", col=c("#aaaaaa")),
    callScatter = list(col=c("#aaaaaa", LOSS="blue", GAIN="red", LOH="purple")),
    smoothScatter = list(pch=".", col=c("#666666")),
    level = list(lty=1L, col=c("black", tcn="purple", c1="blue", c2="red", dh="orange")),
    callLevel = list(lty=1L, col=c("#666666")),
    knownSegment = list(lty=1L, col=c("#aaaaaa"))
  )

  getOptionValue <- function(option, what, track, ...) {
    values <- opts[[option]][[what]]
    value <- values[track]
    if (is.na(value)) value <- values[1L]
    unname(value)
  } # getOptionValue()

  getScatterColor <- function(track, ...) {
    getOptionValue("scatter", "col", track, ...)
  } # getScatterColor()

  getLevelColor <- function(track, ...) {
    getOptionValue("level", "col", track, ...)
  } # getLevelColor()

  getCallScatterColor <- function(track, ...) {
    getOptionValue("callScatter", "col", track, ...)
  } # getCallScatterColor()

  getCallLevelColor <- function(track, ...) {
    getOptionValue("callLevel", "col", track, ...)
  } # getCallLevelColor()

  getCallLevelLty <- function(track, ...) {
    getOptionValue("callLevel", "lty", track, ...)
  } # getCallLevelColor()

  getCIColor <- function(track, ...) {
    getLevelColor(track, ...)
  } # getLevelColor()

  getKnownSegmentColor <- function(track, ...) {
    getOptionValue("knownSegment", "col", track, ...)
  } # getKnownSegmentColor()

  getKnownSegmentLty <- function(track, ...) {
    getOptionValue("knownSegment", "lty", track, ...)
  } # getKnownSegmentColor()

  drawXLabelTicks <- function() {
    if (identical(xlabTicks, "[chr]")) {
      mtext(text=chrTags, side=rep(c(1,3), length.out=length(chrTags)), at=mids, line=0.1, cex=0.7)
    } else if (identical(xlabTicks, "[pos]")) {
      axis(side=1L)
    }
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'fit':

  # Argument 'chromosomes':
  if (!is.null(chromosomes)) {
    disallow <- c("NaN", "Inf")
    chromosomes <- Arguments$getIntegers(chromosomes, range=c(0,Inf), disallow=disallow)
    .stop_if_not(all(is.element(chromosomes, getChromosomes(fit))))
  }

  # Argument 'tracks':
  knownTracks <- c("tcn", "dh", "tcn,c1,c2", "c1,c2", "c1", "c2",
                   "betaN", "betaT", "betaTN")
  defaultTracks <- knownTracks[1:3]
  if (is.null(tracks)) {
    tracks <- defaultTracks
  } else {
    tracks <- match.arg(tracks, choices=knownTracks, several.ok=TRUE)
    tracks <- unique(tracks)
  }

  # Argument 'scatter':
  if (!is.null(scatter)) {
    scatter <- Arguments$getCharacter(scatter)
    if (scatter == "*") {
      scatter <- tracks
    } else {
      scatterT <- strsplit(scatter, split=",", fixed=TRUE)
      tracksT <- strsplit(tracks, split=",", fixed=TRUE)
      .stop_if_not(all(is.element(scatterT, tracksT)))
      # Not needed anymore
      scatterT <- tracksT <- NULL
    }
  }

  # Argument 'calls':
  if (!is.null(calls)) {
    calls <- sapply(calls, FUN=Arguments$getRegularExpression)
  }

  # Argument 'callLoci':
  callLoci <- Arguments$getLogical(callLoci)

  # Argument 'callThresholds':
  callThresholds <- Arguments$getLogical(callThresholds)

  # Argument 'boundaries':
  boundaries <- Arguments$getLogical(boundaries)

  # Argument 'knownSegments':
  knownSegments <- Arguments$getLogical(knownSegments)

  # Argument 'add':
  add <- Arguments$getLogical(add)

  # Argument 'Clim' & 'Blim':
  if (!add) {
    Clim <- Arguments$getNumerics(Clim, length=c(2L,2L),
                                        disallow=c("Inf", "NA", "NaN"))
    Blim <- Arguments$getNumerics(Blim, length=c(2L,2L),
                                        disallow=c("Inf", "NA", "NaN"))
  }

  # Argument 'xScale':
  xScale <- Arguments$getNumeric(xScale, range=c(0,Inf))

  # Argument 'xlabTicks':
  if (!is.null(xlabTicks)) {
    xlabTicks <- Arguments$getCharacter(xlabTicks)
  }

  # Argument 'subset':
  if (!is.null(subset)) {
    subset <- Arguments$getDouble(subset)
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    pushState(verbose)
    on.exit(popState(verbose))
  }


  verbose && enter(verbose, "Plotting PSCN tracks")

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Subset by chromosomes
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!is.null(chromosomes)) {
    verbose && enter(verbose, "Plotting a subset of the chromosomes")
    fit <- extractChromosomes(fit, chromosomes=chromosomes, verbose=verbose)
    verbose && exit(verbose)
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Tile chromosomes
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  fit <- tileChromosomes(fit, verbose=verbose)
  verbose && str(verbose, fit)

  # Extract the input data
  data <- getLocusData(fit)
  if (is.null(data)) {
    throw("Cannot plot segmentation results. No input data available.")
  }

  # Extract the segmentation
  segs <- as.data.frame(fit)

  # Identify available calls
  callData <- NULL
  if (!is.null(calls) || callThresholds) {
    verbose && enter(verbose, "Identifying calls")

    pattern <- "Call$"
    allCallColumns <- grep(pattern, colnames(segs), value=TRUE)
    allCallLabels <- toupper(gsub(pattern, "", allCallColumns))
    verbose && cat(verbose, "Call columns:")
    verbose && print(verbose, allCallColumns)

    if (!is.null(calls)) {
      callColumns <- allCallColumns
      if (length(callColumns) > 0L) {
        keep <- sapply(calls, FUN=function(pattern) {
          (regexpr(pattern, callColumns) != -1L)
        })
        if (is.matrix(keep)) {
          keep <- rowAnys(keep)
        }
        callColumns <- callColumns[keep]
        callLabels <- allCallLabels[keep]

        # Annotate individual loci by calls?
        if (callLoci) {
          callData <- extractCallsByLocus(fit, verbose=less(verbose,5))
        }
      }
      verbose && cat(verbose, "Call to be annotated:")
      verbose && print(verbose, callColumns)
    }

    verbose && exit(verbose)
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Subset of the loci?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!is.null(subset)) {
    # (a) Set and unset the random seed
    if (!is.null(seed)) {
      randomSeed("set", seed=seed, kind="L'Ecuyer-CMRG")
      on.exit(randomSeed("reset"), add=TRUE)
      verbose && printf(verbose, "Random seed temporarily set (seed=%d)\n", seed)
    }

    # (b) Subset
    n <- nrow(data)
    keep <- sample(n, size=subset*n)
    data <- data[keep,]
    if (!is.null(callData)) {
      callData <- callData[keep,]
    }
  }

  # To please R CMD check
  CT <- rho <- muN <- betaT <- betaN <- betaTN <- rho <- NULL
  rm(list=c("CT", "rho", "muN", "betaT", "betaN", "betaTN"))
  attachLocally(data)
  # Calculate (C1,C2)
  C1 <- 1/2*(1-rho)*CT
  C2 <- CT - C1

  # BACKWARD COMPATIBILITY:
  # If 'rho' is not available, recalculate it from tumor BAFs.
  # NOTE: This should throw an error in the future. /HB 2013-10-25
  if (is.null(data$rho)) {
    isSnp <- (!is.na(betaTN) & !is.na(muN))
    isHet <- isSnp & (muN == 1/2)
    rho <- rep(NA_real_, times=nbrOfLoci)
    rho[isHet] <- 2*abs(betaTN[isHet]-1/2)
    warning(sprintf("Locus-level DH signals ('rho') were not available in the %s object and therefore recalculated from the TumorBoost-normalized tumor BAFs ('betaTN').", class(fit)[1L]))
  }

  x <- xScale * x
  vs <- xScale * fit$chromosomeStats[,1:2,drop=FALSE]
  mids <- (vs[,1]+vs[,2])/2

  nbrOfLoci <- length(x)
  chrTags <- sprintf("Chr%02d", chromosomes)

  if (subplots) {
    subplots(length(tracks), ncol=1L)
    par(oma=oma, mar=mar)
  }

  pchT <- if (!is.null(scatter)) { pch } else { NA }
  xlim <- range(x, na.rm=TRUE)
  xlab <- "Genomic position"

  # Graphical handle
  gh <- list(fit=fit)
  gh$xScale <- xScale
  gh$xlim <- xlim
  gh$xlab <- xlab
  if (!is.null(callData)) {
    gh$callsByLocus <- callData
  }

  for (tt in seq_along(tracks)) {
    track <- tracks[tt]
    verbose && enter(verbose, sprintf("Track #%d ('%s') of %d",
                                             tt, track, length(tracks)))

    # Get graphical style parameters.
    tracksT <- unlist(strsplit(track, split=",", fixed=TRUE))
    colS <- sapply(tracksT, FUN=getScatterColor)
    colL <- sapply(tracksT, FUN=getLevelColor)
    colC <- sapply(tracksT, FUN=getCIColor)

    # Color scatter plot according to calls?
    if (!is.null(calls) && callLoci && length(callColumns) > 0L) {
      colsT <- rep(colS[1L], times=nrow(callData))
      for (cc in seq_along(callColumns)) {
        callColumn <- callColumns[cc]
        callLabel <- callLabels[cc]
        verbose && enter(verbose, sprintf("Call #%d ('%s') of %d",
                                      cc, callLabel, length(callColumns)))

        verbose && cat(verbose, "Column: ", callColumn)

        skip <- TRUE
        if (regexpr("tcn", track) != -1L) {
          skip <- !is.element(callLabel, c("LOSS", "NTCN", "GAIN", "LOH"))
        } else if (track == "dh") {
          skip <- !is.element(callLabel, c("AB", "LOH"))
        }
        if (skip) {
          verbose && exit(verbose)
          next
        }

        callsCC <- callData[[callColumn]]
        idxs <- which(callsCC)
        # Nothing to do?
        if (length(idxs) == 0L) {
          verbose && exit(verbose)
          next
        }

        callCol <- getCallScatterColor(callLabel)

        colsT[idxs] <- callCol
      } # for (cc in ...)

      colS <- colsT
    } # if (!is.null(calls))


    # Assign graphical-handle parameters
    gh$track <- track
    gh$subtracks <- tracksT
    gh$scatter <- list(col=colS, pch=pchT)
    gh$level <- list(col=colL)
    gh$cis <- list(col=colC)


    if (track == "tcn") {
      plot(NA, xlim=xlim, ylim=Clim, xlab=xlab, ylab="TCN", axes=FALSE)
      if (!is.null(onBegin)) attachGH(onBegin(gh=gh))
      if (!is.na(pchT)) {
        points(x, CT, pch=pchT, col=colS)
      }
      drawConfidenceBands(fit, what="tcn", quantiles=quantiles, col=colC["tcn"], xScale=xScale)
      drawLevels(fit, what="tcn", col=colL, xScale=xScale)
    }

    if (is.element(track, c("tcn,c1,c2", "c1,c2", "c1", "c2"))) {
      tracksT <- unlist(strsplit(track, split=",", fixed=TRUE))
      plot(NA, xlim=xlim, ylim=Clim, xlab=xlab, ylab="C1, C2, TCN", axes=FALSE)
      if (!is.null(onBegin)) attachGH(onBegin(gh=gh))

      # Draw scatter for TCN or C1 and C2.
      if (!is.na(pchT)) {
        if (is.element("tcn", tracksT)) {
          points(x, CT, pch=pchT, col=colS)
        } else {
          if (is.element("c1", tracksT)) {
            points(x, C1, pch=pchT, col=colS)
          }
          if (is.element("c2", tracksT)) {
            points(x, C2, pch=pchT, col=colS)
          }
        }
      }

      # Draw confidence bands for TCN, C1, C2.
      if (is.element("tcn", tracksT)) {
        drawConfidenceBands(fit, what="tcn", quantiles=quantiles, col=colC["tcn"], xScale=xScale)
      }
      if (is.element("c2", tracksT)) {
        drawConfidenceBands(fit, what="c2", quantiles=quantiles, col=colC["c2"], xScale=xScale)
      }
      if (is.element("c1", tracksT)) {
        drawConfidenceBands(fit, what="c1", quantiles=quantiles, col=colC["c1"], xScale=xScale)
      }

      # Draw segment means for TCN, C1, C2.
      if (is.element("tcn", tracksT)) {
        drawLevels(fit, what="tcn", col=colL["tcn"], xScale=xScale)
      }
      if (is.element("c2", tracksT)) {
        drawLevels(fit, what="c2", col=colL["c2"], xScale=xScale)
      }
      if (is.element("tcn", tracksT)) {
        # In case C2 overlaps with TCN
        drawLevels(fit, what="tcn", col=colL["tcn"], lty="22", xScale=xScale)
      }
      # In case C1 overlaps with C2
      if (is.element("c1", tracksT)) {
        drawLevels(fit, what="c1", col=colL["c1"], xScale=xScale)
        if (is.element("c2", tracksT)) {
          drawLevels(fit, what="c2", col=colL["c2"], lty="22", xScale=xScale)
        }
        if (is.element("tcn", tracksT)) {
          drawLevels(fit, what="tcn", col=colL["tcn"], lty="22", xScale=xScale)
        }
      }
    }

    if (track == "betaN") {
      plot(NA, xlim=xlim, ylim=Blim, xlab=xlab, ylab="BAF_N", axes=FALSE)
      if (!is.null(onBegin)) attachGH(onBegin(gh=gh))
      if (!is.na(pchT)) {
        points(x, betaN, pch=pchT, col=colS)
      }
    }

    if (track == "betaT") {
      plot(NA, xlim=xlim, ylim=Blim, xlab=xlab, ylab="BAF_T", axes=FALSE)
      if (!is.null(onBegin)) attachGH(onBegin(gh=gh))
      if (!is.na(pchT)) {
        points(x, betaT, pch=pchT, col=colS)
      }
    }

    if (track == "betaTN") {
      plot(NA, xlim=xlim, ylim=Blim, xlab=xlab, ylab="BAF_TN", axes=FALSE)
      if (!is.null(onBegin)) attachGH(onBegin(gh=gh))
      if (!is.na(pchT)) {
        points(x, betaTN, pch=pchT, col=colS)
      }
    }

    if (track == "dh") {
      plot(NA, xlim=xlim, ylim=Blim, xlab=xlab, ylab="DH", axes=FALSE)
      if (!is.null(onBegin)) attachGH(onBegin(gh=gh))
      if (!is.na(pchT)) {
        points(x, rho, pch=pchT, col=colS)
      }
      drawConfidenceBands(fit, what="dh", quantiles=quantiles, col=colC["dh"], xScale=xScale)
      drawLevels(fit, what="dh", col=colL["dh"], xScale=xScale)
    }

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # For each panel of tracks, annotate segments with calls?
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (!is.null(calls) && !callLoci && length(callColumns) > 0L) {
      for (cc in seq_along(callColumns)) {
        callColumn <- callColumns[cc]
        callLabel <- callLabels[cc]
        verbose && enter(verbose, sprintf("Call #%d ('%s') of %d",
                                      cc, callLabel, length(callColumns)))

        verbose && cat(verbose, "Column: ", callColumn)

        segsT <- segs[,c("dhStart", "dhEnd", callColumn)]
        isCalled <- which(segsT[[callColumn]])
        segsT <- segsT[isCalled,1:2,drop=FALSE]
        verbose && printf(verbose, "Number of segments called %s: %d\n",
                                                  callLabel, nrow(segsT))
        segsT <- xScale * segsT

        verbose && str(verbose, segsT)

        side <- 2*((cc+1) %% 2) + 1
        # For each segment called...
        for (ss in seq_len(nrow(segsT))) {
          x0 <- segsT[ss,1,drop=TRUE]
          x1 <- segsT[ss,2,drop=TRUE]
          abline(v=c(x0,x1), lty=3, col="gray")
          xMid <- (x0+x1)/2
          mtext(side=side, at=xMid, line=-1, cex=0.7, col="#666666", callLabel)
        } # for (ss in ...)
        verbose && exit(verbose)
      } # for (cc in ...)
    } # if (!is.null(calls))


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # For each panel of tracks, annotate with call thresholds?
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (callThresholds) {
      # Add call parameter estimates, e.g. deltaAB
      colCL <- sapply(tracksT, FUN=getCallLevelColor)
      ltyCL <- sapply(tracksT, FUN=getCallLevelLty)

      trackT <- track
      for (cc in seq_along(allCallColumns)) {
        callColumn <- allCallColumns[cc]
        callLabel <- allCallLabels[cc]

        h <- NULL
        if (callLabel == "AB") {
          if (track == "dh") {
            h <- fit$params$deltaAB
            label <- expression(Delta[AB])
          }
        } else if (callLabel == "LOH") {
          if (regexpr("c1", track) != -1L) {
            h <- fit$params$deltaLowC1
            label <- expression(Delta[LOH])
            trackT <- "c1"
          }
        } else if (callLabel == "NTCN") {
          if (track == "tcn") {
            h <- fit$params$ntcnRange
            label <- c(expression(Delta[-NTCN]), expression(Delta[+NTCN]))
          }
        }

        if (!is.null(h)) {
          abline(h=h, lty=ltyCL[trackT], lwd=2, col=colCL[trackT])
          for (ss in 1:2) {
            side <- c(2,4)[ss]
            adj <- c(1.2,-0.2)[ss]
            mtext(side=side, at=h, label, adj=adj, las=2, xpd=TRUE)
          }
        }
      } # for (cc in ...)
    } # if (callThresholds)

    drawXLabelTicks()
    if (boundaries) {
      abline(v=vs, lty=1, lwd=2)
    }
    if (knownSegments) {
      colT <- getKnownSegmentColor()
      ltyT <- getKnownSegmentLty()
      drawKnownSegments(fit, col=colT, lty=ltyT)
    }
    axis(side=2); box()
    if (!is.null(onEnd)) onEnd(gh=gh)

    verbose && exit(verbose)
  } # for (tt ...)

  verbose && exit(verbose)

  invisible(gh)
}, private=TRUE) # plotTracksManyChromosomes()
