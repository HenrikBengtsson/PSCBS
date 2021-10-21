###########################################################################/**
# @set "class=PairedPSCBS"
# @RdocMethod plotTracks
# @aliasmethod plotTracks1
# @aliasmethod plotTracks2
# @aliasmethod plotTracksManyChromosomes
# @alias plot
#
# @title "Plots parental specific copy numbers along the genome"
#
# \description{
#  @get "title" for one or more chromosomes.
#  It is possible to specify what type of tracks to plot.
#  Each type of track is plotted in its own panel.
# }
#
# @synopsis
#
# \arguments{
#   \item{x}{A result object returned by @see "segmentByPairedPSCBS".}
#   \item{tracks}{A @character @vector specifying what types of tracks to plot.}
#   \item{scatter}{A @character @vector specifying which of the tracks should
#     have scatter plot.}
#   \item{calls}{A @character @vector of regular expression identifying
#     call labels to be highlighted in the panels.}
#   \item{pch}{The type of the scatter points, if any.}
#   \item{col}{The color of the scatter points, if any.}
#   \item{cex}{The size of the scatter points, if any.}
#   \item{changepoints}{If @TRUE, changepoints are drawn as vertical lines.}
#   \item{grid}{If @TRUE, horizontal lines are displayed.}
#   \item{quantiles}{A @numeric @vector in [0,1] specifying the quantiles
#      of the confidence bands to be drawn, if any.}
#   \item{xlim}{(Optional) The genomic range to plot.}
#   \item{Clim}{The range of copy numbers.}
#   \item{Blim}{The range of allele B fractions (BAFs) and
#     decrease of heterozygosity (DHs).}
#   \item{xScale}{The scale factor used for genomic positions.}
#   \item{...}{Not used.}
#   \item{add}{If @TRUE, the panels plotted are added to the existing plot,
#     otherwise a new plot is created.}
#   \item{subplots}{If @TRUE, then subplots are automatically setup.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   Returns nothing.
# }
#
# @author "HB"
#
# @keyword IO
# @keyword internal
#*/###########################################################################
setMethodS3("plotTracks1", "PairedPSCBS", function(x, tracks=c("tcn", "dh", "tcn,c1,c2", "tcn,c1", "tcn,c2", "c1,c2", "betaN", "betaT", "betaTN")[1:3], scatter="*", calls=".*", pch=".", col=NULL, cex=1, changepoints=FALSE, grid=FALSE, quantiles=c(0.05,0.95), xlim=NULL, Clim=c(0,3*ploidy(x)), Blim=c(0,1), xScale=1e-6, ..., add=FALSE, subplots=!add && (length(tracks) > 1), verbose=FALSE) {

  # To please R CMD check
  fit <- x

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'add':
  add <- Arguments$getLogical(add)

  # Argument 'Clim' & 'Blim':
  if (!add) {
    Clim <- Arguments$getNumerics(Clim, length=c(2L,2L),
                                        disallow=c("Inf", "NA", "NaN"))
    Blim <- Arguments$getNumerics(Blim, length=c(2L,2L),
                                        disallow=c("Inf", "NA", "NaN"))
  }

  # Argument 'fit':
  if (nbrOfChromosomes(fit) >= 1L) {
    return(plotTracksManyChromosomes(fit, tracks=tracks, scatter=scatter, calls=calls, pch=pch, quantiles=quantiles, Clim=Clim, Blim=Blim, xScale=xScale, ..., add=add, subplots=subplots, verbose=verbose))
  }

  # Argument 'tracks':
  knownTracks <- c("tcn", "dh", "tcn,c1,c2", "tcn,c1", "tcn,c2", "c1,c2", "betaN", "betaT", "betaTN")
  tracks <- match.arg(tracks, choices=knownTracks, several.ok=TRUE)
  tracks <- unique(tracks)

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

  # Argument 'changepoints':
  changepoints <- Arguments$getLogical(changepoints)

  # Argument 'grid':
  grid <- Arguments$getLogical(grid)

  # Argument 'xlim':
  if (!is.null(xlim)) {
    xlim <- Arguments$getNumerics(xlim, length=c(2,2))
  }

  # Argument 'xScale':
  xScale <- Arguments$getNumeric(xScale, range=c(0,Inf))

  # Argument 'subplots':
  subplots <- Arguments$getLogical(subplots)

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    pushState(verbose)
    on.exit(popState(verbose))
  }


  verbose && enter(verbose, "Plotting PSCN tracks")

  # Extract the input data
  data <- getLocusData(fit)
  if (is.null(data)) {
    stop("Cannot plot segmentation results. No input data available.")
  }

  chromosomes <- getChromosomes(fit)
  chromosome <- chromosomes[1]
  x <- data$x
  CT <- data$CT
  rho <- data$rho
  betaT <- data$betaT
  betaN <- data$betaN
  betaTN <- data$betaTN
  muN <- data$muN
  rho <- data$rho
  hasDH <- !is.null(rho)
  if (hasDH) {
    isHet <- !is.na(rho)
    isSnp <- isHet
  } else {
    isSnp <- (!is.na(betaTN) & !is.na(muN))
    isHet <- isSnp & (muN == 1/2)
  }
  nbrOfLoci <- length(x)

  # BACKWARD COMPATIBILITY:
  # If 'rho' is not available, recalculate it from tumor BAFs.
  # NOTE: This should throw an error in the future. /HB 2013-10-25
  if (!hasDH) {
    rho <- rep(NA_real_, times=nbrOfLoci)
    rho[isHet] <- 2*abs(betaTN[isHet]-1/2)
    warning(sprintf("Locus-level DH signals ('rho') were not available in the %s object and therefore recalculated from the TumorBoost-normalized tumor BAFs ('betaTN').", class(fit)[1L]))
  }

  # Extract the segmentation
  segs <- as.data.frame(fit)

  # Identify available calls
  if (!is.null(calls)) {
    verbose && enter(verbose, "Identifying calls")

    pattern <- "Call$"
    callColumns <- grep(pattern, colnames(segs), value=TRUE)
    if (length(callColumns) > 0) {
      keep <- sapply(calls, FUN=function(pattern) {
        (regexpr(pattern, callColumns) != -1)
      })
      if (is.matrix(keep)) {
        keep <- rowAnys(keep, useNames=FALSE)
      }
      callColumns <- callColumns[keep]
      callLabels <- gsub(pattern, "", callColumns)
      callLabels <- toupper(callLabels)
    }
    verbose && cat(verbose, "Call columns:")
    verbose && print(verbose, callColumns)

    verbose && exit(verbose)
  } else {
    callColumns <- NULL
  }

  if (chromosome != 0) {
    chrTag <- sprintf("Chr%02d", chromosome)
  } else {
    chrTag <- ""
  }

  if (xScale != 1) {
    x <- xScale * x
    if (!is.null(xlim)) {
      xlim <- xScale * xlim
    }
  }

  if (subplots) {
    subplots(length(tracks), ncol=1)
    par(mar=c(1,4,1,2)+1)
  }

  # Color loci by heterozygous vs homozygous
  if (hasDH) {
    colMu <- c("gray", "black")[!is.na(rho) + 1]
  } else {
    colMu <- c("gray", "black")[(muN == 1/2) + 1]
  }

  for (tt in seq_along(tracks)) {
    track <- tracks[tt]
    verbose && enter(verbose, sprintf("Track #%d ('%s') of %d",
                                             tt, track, length(tracks)))

    if (!is.null(scatter)) {
      pchT <- pch
      colT <- col
    } else {
      pchT <- NA
      colT <- NA
    }

    if (track == "tcn") {
      colT <- ifelse(is.null(colT), "black", colT)
      if (add) {
        points(x, CT, pch=pchT, col=colT, cex=cex)
      } else {
        plot(x, CT, pch=pchT, col=colT, cex=cex, xlim=xlim, ylim=Clim, ylab="TCN")
        stext(side=3, pos=1, chrTag)
        if (grid) {
          abline(h=seq(from=0, to=Clim[2], by=2), lty=3, col="gray")
          abline(h=0, lty=1, col="black")
        }
        drawConfidenceBands(fit, what="tcn", quantiles=quantiles, col="purple", xScale=xScale)
        drawLevels(fit, what="tcn", col="purple", xScale=xScale)
      }
    }

    if (is.element(track, c("tcn,c1,c2", "tcn,c1", "tcn,c2", "c1,c2"))) {
      colT <- ifelse(is.null(colT), "black", colT)
      subtracks <- strsplit(track, split=",", fixed=TRUE)[[1]]
      ylab <- paste(toupper(subtracks), collapse=", ")
      if (add) {
        points(x, CT, pch=pchT, cex=cex, col=colT)
      } else {
        plot(x, CT, pch=pchT, cex=cex, col=colT, xlim=xlim, ylim=Clim, ylab=ylab)
        stext(side=3, pos=1, chrTag)
        if (grid) {
          abline(h=seq(from=0, to=Clim[2], by=2), lty=3, col="gray")
          abline(h=0, lty=1, col="black")
        }
        if (is.element("tcn", subtracks)) {
          drawConfidenceBands(fit, what="tcn", quantiles=quantiles, col="purple", xScale=xScale)
        }
        if (is.element("c2", subtracks)) {
          drawConfidenceBands(fit, what="c2", quantiles=quantiles, col="red", xScale=xScale)
        }
        if (is.element("c1", subtracks)) {
          drawConfidenceBands(fit, what="c1", quantiles=quantiles, col="blue", xScale=xScale)
        }
        if (is.element("tcn", subtracks)) {
          drawLevels(fit, what="tcn", col="purple", xScale=xScale)
        }
        if (is.element("c2", subtracks)) {
          drawLevels(fit, what="c2", col="red", xScale=xScale)
          # In case C2 overlaps with TCN
          if (is.element("tcn", subtracks)) {
            drawLevels(fit, what="tcn", col="purple", lty="22", xScale=xScale)
          }
        }
        if (is.element("c1", subtracks)) {
          drawLevels(fit, what="c1", col="blue", xScale=xScale)
          # In case C1 overlaps with C1
          if (is.element("c2", subtracks)) {
            drawLevels(fit, what="c2", col="red", lty="22", xScale=xScale)
            if (is.element("tcn", subtracks)) {
              drawLevels(fit, what="tcn", col="purple", lty="22", xScale=xScale)
            }
          }
        }
      }
    }

    if (track == "betaN") {
      colT <- ifelse(is.null(colT), colMu, colT)
      if (add) {
        points(x, betaN, pch=pchT, cex=cex, col="black")
      } else {
        plot(x, betaN, pch=pchT, cex=cex, col=colT, xlim=xlim, ylim=Blim, ylab=expression(BAF[N]))
        stext(side=3, pos=1, chrTag)
      }
    }

    if (track == "betaT") {
      colT <- ifelse(is.null(colT), colMu, colT)
      if (add) {
        points(x, betaT, pch=pchT, cex=cex, col="black")
      } else {
        plot(x, betaT, pch=pchT, cex=cex, col=colT, xlim=xlim, ylim=Blim, ylab=expression(BAF[T]))
        stext(side=3, pos=1, chrTag)
      }
    }

    if (track == "betaTN") {
      colT <- ifelse(is.null(colT), colMu, colT)
      if (add) {
        points(x, betaTN, pch=pchT, cex=cex, col="black")
      } else {
        plot(x, betaTN, pch=pchT, cex=cex, col=colT, xlim=xlim, ylim=Blim, ylab=expression(BAF[T]^"*"))
        stext(side=3, pos=1, chrTag)
      }
    }

    if (track == "dh") {
      colT <- ifelse(is.null(colT), colMu[isHet], colT)
      if (add) {
        points(x, rho, pch=pchT, cex=cex, col="black")
      } else {
        plot(x, rho, pch=pchT, cex=cex, col=colT, xlim=xlim, ylim=Blim, ylab="DH")
        stext(side=3, pos=1, chrTag)
        drawConfidenceBands(fit, what="dh", quantiles=quantiles, col="orange", xScale=xScale)
        drawLevels(fit, what="dh", col="orange", xScale=xScale)
      }
    }

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Draw change points?
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (changepoints) {
      drawChangePoints(fit, col="#666666", xScale=xScale)
    }


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # For each panel of tracks, annotate with calls?
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (length(callColumns) > 0) {
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

      # Add call parameter estimates, e.g. deltaAB
      for (cc in seq_along(callColumns)) {
        callColumn <- callColumns[cc]
        callLabel <- callLabels[cc]
        h <- NULL
        if (callLabel == "AB") {
          if (track == "dh") {
            h <- fit$params$deltaAB
            label <- expression(Delta[AB])
            colT <- "orange"
          }
        } else if (callLabel == "LOH") {
          if (regexpr("c1", track) != -1L) {
            h <- fit$params$deltaLowC1
            label <- expression(Delta[LOH])
            colT <- "blue"
          }
        } else if (callLabel == "NTCN") {
          if (track == "tcn") {
            h <- fit$params$ntcnRange
            label <- c(expression(Delta[-NTCN]), expression(Delta[+NTCN]))
            colT <- "purple"
          }
        }

        if (!is.null(h)) {
          abline(h=h, lty=4, lwd=2, col=colT)
          for (ss in 1:2) {
            side <- c(2,4)[ss]
            adj <- c(1.2,-0.2)[ss]
            mtext(side=side, at=h, label, adj=adj, las=2, xpd=TRUE)
          }
        }
      } # for (cc in ...)
    } # if (length(callColumns) > 0)

    verbose && exit(verbose)
  } # for (tt ...)

  verbose && exit(verbose)

  invisible()
}, private=TRUE) # plotTracks1()


setMethodS3("plotTracks", "PairedPSCBS", function(fit, ...) {
  plotTracksManyChromosomes(fit, ...)
})

setMethodS3("plot", "PairedPSCBS", function(x, ...) {
  plotTracks(x, ...)
}, private=TRUE)


setMethodS3("drawLevels", "PairedPSCBS", function(fit, what=c("tcn", "betaTN", "dh", "c1", "c2"), lend=1L, xScale=1e-6, ...) {
  # WORKAROUND: If Hmisc is loaded after R.utils, it provides a buggy
  # capitalize() that overrides the one we want to use. Until PSCBS
  # gets a namespace, we do the following workaround. /HB 2011-07-14
  capitalize <- R.utils::capitalize

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'what':
  what <- match.arg(what)

  # Argument 'xScale':
  xScale <- Arguments$getNumeric(xScale, range=c(0,Inf))


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Tile chromosomes
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  fitT <- tileChromosomes(fit)

  # Get segmentation results
  segs <- as.data.frame(fitT)

  if (what == "betaTN") {
    whatT <- "dh"
  } else {
    whatT <- what
  }

  # Extract subset of segments
  fields <- c("start", "end")
  fields <- sprintf("%s%s", ifelse(what == "tcn", what, "dh"), capitalize(fields))
  fields <- c(fields, sprintf("%sMean", whatT))
  segsT <- segs[,fields, drop=FALSE]
  segsT <- unique(segsT)

  if (what == "betaTN") {
    dh <- segsT[,"dhMean"]
    bafU <- (1 + dh)/2
    bafL <- (1 - dh)/2
    segsT[,3] <- bafU
    segsT[,4] <- bafL
  }

  # Reuse drawLevels() for the DNAcopy class
  for (cc in seq(from=3, to=ncol(segsT))) {
    segsTT <- segsT[,c(1:2,cc)]
    colnames(segsTT) <- c("loc.start", "loc.end", "seg.mean")
    dummy <- list(output=segsTT)
    class(dummy) <- "DNAcopy"
    drawLevels(dummy, lend=lend, xScale=xScale, ...)
  } # for (cc ...)
}, private=TRUE)


setMethodS3("drawConfidenceBands", "PairedPSCBS", function(fit, what=c("tcn", "dh", "c1", "c2"), quantiles=c(0.05,0.95), col=col, alpha=0.4, xScale=1e-6, ...) {
  # WORKAROUND: If Hmisc is loaded after R.utils, it provides a buggy
  # capitalize() that overrides the one we want to use. Until PSCBS
  # gets a namespace, we do the following workaround. /HB 2011-07-14
  capitalize <- R.utils::capitalize

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'what':
  what <- match.arg(what)

  # Argument 'quantiles':
  if (!is.null(quantiles)) {
    quantiles <- Arguments$getNumerics(quantiles, range=c(0,1), length=c(2,2))
  }

  # Argument 'xScale':
  xScale <- Arguments$getNumeric(xScale, range=c(0,Inf))


  # Nothing todo?
  if (is.null(quantiles)) {
    return()
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Tile chromosomes
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  fitT <- tileChromosomes(fit)

  # Get segmentation results
  segs <- as.data.frame(fitT)

  # Extract subset of segments
  fields <- c("start", "end")
  fields <- sprintf("%s%s", ifelse(what == "tcn", what, "dh"), capitalize(fields))

  tags <- sprintf("%g%%", 100*quantiles)
  qFields <- sprintf("%s_%s", what, tags)

  # Nothing todo?
  if (!all(is.element(qFields, colnames(segs)))) {
    return()
  }

  fields <- c(fields, qFields)

  segsT <- segs[,fields, drop=FALSE]
  segsT <- unique(segsT)

  # Rescale x-axis
  segsT[,1:2] <- xScale * segsT[,1:2]

  colQ <- col2rgb(col, alpha=TRUE)
  colQ["alpha",] <- alpha*colQ["alpha",]
  colQ <- rgb(red=colQ["red",], green=colQ["green",], blue=colQ["blue",], alpha=colQ["alpha",], maxColorValue=255)

  for (kk in seq_len(nrow(segsT))) {
    rect(xleft=segsT[kk,1], xright=segsT[kk,2], ybottom=segsT[kk,3], ytop=segsT[kk,4], col=colQ, border=FALSE)
  }
}, private=TRUE)



setMethodS3("plotC1C2", "PairedPSCBS", function(fit, ..., xlab=expression(C[1]), ylab=expression(C[2]), Clim=c(0,2*ploidy(fit))) {
  # Argument 'Clim':
  Clim <- Arguments$getNumerics(Clim, length=c(2L,2L),
                                        disallow=c("Inf", "NA", "NaN"))

  plot(NA, xlim=Clim, ylim=Clim, xlab=xlab, ylab=ylab)
  abline(a=0, b=1, lty=3)
  pointsC1C2(fit, ...)
}, private=TRUE)


setMethodS3("pointsC1C2", "PairedPSCBS", function(fit, cex=NULL, col="#00000033", ...) {
  data <- extractC1C2(fit)
  X <- data[,1:2,drop=FALSE]
  n <- data[,4,drop=TRUE]
  n <- sqrt(n)
  w <- n / sum(n, na.rm=TRUE)

  if (is.null(cex)) {
    cex <- w
    cex <- cex / mean(cex, na.rm=TRUE)
    cex <- cex + 1/2
  }

  points(X, cex=cex, col=col, ...)
}, private=TRUE)


setMethodS3("linesC1C2", "PairedPSCBS", function(fit, ...) {
  drawChangePointsC1C2(fit, ...)
}, private=TRUE)


setMethodS3("drawChangePointsC1C2", "PairedPSCBS", function(fit, col="#00000033", labels=FALSE, lcol="#333333", cex=0.7, adj=c(+1.5,+0.5), ...) {
  xy <- extractMinorMajorCNs(fit, splitters=TRUE, addGaps=TRUE)
  xy <- xy[,1:2,drop=FALSE]
  res <- lines(xy, col=col, ...)

  if (labels) {
    n <- nrow(xy)
    dxy <- (xy[-1,] - xy[-n,]) / 2
    xyMids <- xy[-n,] + dxy
    labels <- rownames(xy)
    labels <- sprintf("%s-%s", labels[-n], labels[-1])
    text(xyMids, labels, cex=cex, col=lcol, adj=adj, ...)
  }

  invisible(res)
}, private=TRUE)




setMethodS3("plotDeltaC1C2", "PairedPSCBS", function(fit, ..., xlab=expression(Delta*C[1]), ylab=expression(Delta*C[2]), Clim=c(-1,1)*ploidy(fit)) {
  # Argument 'Clim':
  Clim <- Arguments$getNumerics(Clim, length=c(2L,2L),
                                        disallow=c("Inf", "NA", "NaN"))

  plot(NA, xlim=Clim, ylim=Clim, xlab=xlab, ylab=ylab)
  abline(h=0, lty=3)
  abline(v=0, lty=3)
  pointsDeltaC1C2(fit, ...)
}, private=TRUE)


setMethodS3("pointsDeltaC1C2", "PairedPSCBS", function(fit, ...) {
  data <- extractDeltaC1C2(fit)
  X <- data[,1:2,drop=FALSE]
  points(X, ...)
}, private=TRUE)


setMethodS3("linesDeltaC1C2", "PairedPSCBS", function(fit, ...) {
  xy <- extractDeltaC1C2(fit)
  xy <- xy[,1:2,drop=FALSE]
  lines(xy, ...)
}, private=TRUE)



setMethodS3("arrowsC1C2", "PairedPSCBS", function(fit, length=0.05, ...) {
  xy <- extractMinorMajorCNs(fit, splitters=TRUE, addGaps=TRUE)
  xy <- xy[,1:2,drop=FALSE]
  x <- xy[,1,drop=TRUE]
  y <- xy[,2,drop=TRUE]
  s <- seq_len(length(x)-1)
  arrows(x0=x[s],y0=y[s], x1=x[s+1],y1=y[s+1], code=2, length=length, ...)
}, private=TRUE)


setMethodS3("arrowsDeltaC1C2", "PairedPSCBS", function(fit, length=0.05, ...) {
  xy <- extractDeltaC1C2(fit)
  xy <- xy[,1:2,drop=FALSE]
  x <- xy[,1,drop=TRUE]
  y <- xy[,2,drop=TRUE]
  s <- seq_len(length(x)-1)
  arrows(x0=x[s],y0=y[s], x1=x[s+1],y1=y[s+1], code=2, length=length, ...)
}, private=TRUE)




setMethodS3("tileChromosomes", "PairedPSCBS", function(fit, chrStarts=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'chrStarts':
  if (!is.null(chrStarts)) {
    chrStarts <- Arguments$getDoubles(chrStarts)
  }

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


  verbose && enter(verbose, "Tile chromosomes")

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
  if (is.null(chrStarts)) {
    xRange <- matrix(0, nrow=length(chromosomes), ncol=2)
    for (kk in seq_along(chromosomes)) {
      chromosome <- chromosomes[kk]
      idxs <- which(data$chromosome == chromosome)
      x <- data$x[idxs]
      r <- range(x, na.rm=TRUE)
      r <- r / 1e6
      r[1] <- floor(r[1])
      r[2] <- ceiling(r[2])
      r <- 1e6 * r
      xRange[kk,] <- r
    } # for (kk ...)

    chrLength <- xRange[,2]
    chrStarts <- c(0, cumsum(chrLength)[-length(chrLength)])
    chrEnds <- chrStarts + chrLength

    # Not needed anymore
    x <- idxs <- NULL
  } # if (is.null(chrStarts))

  verbose && cat(verbose, "Chromosome starts:")
  chromosomeStats <- cbind(start=chrStarts, end=chrEnds, length=chrEnds-chrStarts)
  verbose && print(chromosomeStats)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Offset...
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  segFields <- grep("(Start|End)$", colnames(segs), value=TRUE)
  # Sanity check
  .stop_if_not(length(segFields) > 0)

  for (kk in seq_along(chromosomes)) {
    chromosome <- chromosomes[kk]
    chrTag <- sprintf("Chr%02d", chromosome)
    verbose && enter(verbose, sprintf("Chromosome #%d ('%s') of %d",
                                         kk, chrTag, length(chromosomes)))

    # Get offset for this chromosome
    offset <- chrStarts[kk]
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
  fit$data <- data
  fit$output <- segs
  fit$chromosomeStats <- chromosomeStats
  fit$params$knownSegments <- knownSegments
#  fitT$params$chrOffsets <- chrOffsets

  # Flag object
  attr(fit, "tiledChromosomes") <- TRUE

  verbose && exit(verbose)

  fit
}, private=TRUE) # tileChromosomes()



setMethodS3("drawChangePoints", "PSCBS", function(fit, labels=FALSE, col="#666666", cex=0.5, xScale=1e-6, side=3, line=-1, xpd=TRUE, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Tile chromosomes
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  fitT <- tileChromosomes(fit)
  verbose && str(verbose, fitT)

  segs <- getSegments(fitT, splitters=FALSE)
  xStarts <- segs[,"tcnStart"]
  xEnds <- segs[,"tcnEnd"]

  xs <- sort(unique(c(xStarts, xEnds)))
  abline(v=xScale*xs, lty=1, col=col)

  if (labels) {
    xMids <- xScale * (xEnds + xStarts) / 2
    labels <- rownames(segs)
    mtext(side=side, at=xMids, labels, line=line, cex=cex, col=col, xpd=xpd, ...)
  }
}, protected=TRUE)



setMethodS3("getChromosomeRanges", "PairedPSCBS", function(fit, ...) {
  # To please R CMD check, cf. subset()
  chromosome <- NULL; rm(list="chromosome")

  segs <- getSegments(fit, splitters=FALSE)
  chromosomes <- sort(unique(segs$chromosome))

  # Allocate
  naValue <- NA_real_
  res <- matrix(naValue, nrow=length(chromosomes), ncol=3)
  rownames(res) <- chromosomes
  colnames(res) <- c("start", "end", "length")

  # Get start and end of each chromosome.
  for (ii in seq_len(nrow(res))) {
    chr <- chromosomes[ii]
    segsII <- subset(segs, chromosome == chr)
    res[ii,"start"] <- min(segsII$tcnStart, na.rm=TRUE)
    res[ii,"end"] <- max(segsII$tcnEnd, na.rm=TRUE)
  } # for (ii ...)

  res[,"length"] <- res[,"end"] - res[,"start"] + 1L

  # Sanity check
  .stop_if_not(nrow(res) == length(chromosomes))

  res <- as.data.frame(res)
  res <- cbind(chromosome=chromosomes, res)

  res
}, protected=TRUE) # getChromosomeRanges()


setMethodS3("getChromosomeOffsets", "PairedPSCBS", function(fit, resolution=1e6, ...) {
  # Argument 'resolution':
  if (!is.null(resolution)) {
    resolution <- Arguments$getDouble(resolution, range=c(1,Inf))
  }

  data <- getChromosomeRanges(fit, ...)
  splits <- data[,"start"] + data[,"length"]

  if (!is.null(resolution)) {
    splits <- ceiling(splits / resolution)
    splits <- resolution * splits
  }

  offsets <- c(0L, cumsum(splits))
  names(offsets) <- c(rownames(data), NA)

  offsets
}, protected=TRUE) # getChromosomeOffsets()
