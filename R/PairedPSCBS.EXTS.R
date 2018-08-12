setMethodS3("shiftTCN", "PairedPSCBS", function(fit, shift, update=TRUE, ...) {
  # Argument 'shift':
  shift <- Arguments$getDouble(shift, disallow=c("NA", "NaN", "Inf"))

  data <- getLocusData(fit)
  data$CT <- data$CT + shift
  fit$data <- data
  # Not needed anymore
  data <- NULL

  if (update) {
    fit <- updateMeans(fit, ...)
  }

  fit
}, protected=TRUE)


setMethodS3("bootstrapCIs", "PairedPSCBS", function(fit, ...) {
  # ...
}, private=TRUE)


###########################################################################/**
# @set "class=PairedPSCBS"
# @RdocMethod extractTCNAndDHs
#
# @title "Extract TCN and DH mean levels per segment"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Arguments passed to \code{getSegments()}.}
# }
#
# \value{
#   Returns a @data.frame.
# }
#
# @author "HB"
#
# \seealso{
#   @seemethod "extractMinorMajorCNs".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("extractTCNAndDHs", "PairedPSCBS", function(fit, ...) {
  segs <- getSegments(fit, ...)
  .stop_if_not(!is.null(segs))

  data <- segs[,c("tcnMean", "dhMean", "tcnNbrOfLoci", "dhNbrOfLoci"), drop=FALSE]
  data
}, protected=TRUE)



###########################################################################/**
# @set "class=PairedPSCBS"
# @RdocMethod extractMinorMajorCNs
# @aliasmethod extractC1C2
#
# @title "Extract minor and major copy-number mean levels per segment"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a @data.frame.
# }
#
# @author "HB"
#
# \seealso{
#   @seemethod "extractTCNAndDHs"
#   @seeclass
# }
#*/###########################################################################
setMethodS3("extractMinorMajorCNs", "PairedPSCBS", function(fit, ...) {
  data <- extractTCNAndDHs(fit, ...)

  gamma <- data[,1L]
  rho <- data[,2L]
  C1 <- 1/2*(1-rho)*gamma
  C2 <- gamma - C1

  data[,1L] <- C1
  data[,2L] <- C2
  colnames(data)[1:2] <- c("C1", "C2")

  # Swap (C1,C2)?
  segs <- getSegments(fit, ...)
  flipped <- segs$c1c2Swap
  if (!is.null(flipped)) {
    idxs <- which(flipped)
    if (length(idxs) > 0L) {
      data[idxs,1:2] <- data[idxs,2:1]
    }
  }

  data
}, protected=TRUE)


setMethodS3("extractC1C2", "PairedPSCBS", function(...) {
  extractMinorMajorCNs(...)
}, protected=TRUE)


setMethodS3("extractCNs", "PairedPSCBS", function(fit, splitters=TRUE, ...) {
  data <- extractC1C2(fit, splitters=splitters, ...)
  data[,c("C1", "C2"), drop=FALSE]
})


setMethodS3("extractDeltaC1C2", "PairedPSCBS", function(...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  X <- extractC1C2(..., splitters=TRUE, addGaps=TRUE)

  # (C1,C2)
  C1C2 <- X[,1:2,drop=FALSE]

  # Number of TCN and DH data points
  counts <- X[,3:4, drop=FALSE]

  # Not needed anymore
  X <- NULL

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Calculate (dC1,dC2)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # (dC1, dC2)
  dC1C2 <- matrixStats::colDiffs(C1C2)

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Change-point weights
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Region weights from DH counts
  w <- counts[,2,drop=TRUE]
  w <- sqrt(w)
  w <- w / sum(w, na.rm=TRUE)

  # (a) Smallest of the two flanking (DH) counts
  cpw <- cbind(w[1:(length(w)-1)], w[2:length(w)])
  cpw <- rowMins(cpw, na.rm=TRUE)
  cpw[is.infinite(cpw)] <- NA
  cpw <- sqrt(cpw)
  cpwMin <- cpw / sum(cpw, na.rm=TRUE)

  # (b) Sum of region weights
  cpw <- w[1:(length(w)-1)] + w[2:length(w)]
  cpwAvg <- cpw / sum(cpw, na.rm=TRUE)

  cbind(dC1=dC1C2[,1], dC2=dC1C2[,2], wMin=cpwMin, wAvg=cpwAvg)
}, protected=TRUE)



setMethodS3("postsegmentTCN", "PairedPSCBS", function(fit, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    pushState(verbose)
    on.exit(popState(verbose))
  }

  verbose && enter(verbose, "Post-segmenting TCNs")

  flavor <- fit$params$flavor
  if (!force && regexpr("&", flavor, fixed=TRUE) != -1) {
    verbose && cat(verbose, "Nothing to do. Already postsegmentTCN:ed: ", flavor)
    verbose && exit(verbose)
    return(fit)
  }

  joinSegments <- fit$params$joinSegments
  if (!joinSegments) {
    throw("Postsegmentation of TCNs is only implemented for the case when joinSegments=TRUE: ", joinSegments)
  }


  # Get mean estimators
  estList <- getMeanEstimators(fit, "tcn")
  avgTCN <- estList$tcn


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract the data and segmentation results
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  data <- getLocusData(fit)

  segs <- getSegments(fit)
  keep <- is.finite(segs$chromosome)
  segs <- segs[keep,,drop=FALSE]
  tcnSegRows <- fit$tcnSegRows[keep,,drop=FALSE]
  dhSegRows <- fit$dhSegRows[keep,,drop=FALSE]

  # Sanity check
  .stop_if_not(nrow(dhSegRows) == nrow(tcnSegRows))
  .stop_if_not(all(tcnSegRows[,1] <= tcnSegRows[,2], na.rm=TRUE))
#  .stop_if_not(all(tcnSegRows[-nrow(tcnSegRows),2] < tcnSegRows[-1,1], na.rm=TRUE))
  .stop_if_not(all(dhSegRows[,1] <= dhSegRows[,2], na.rm=TRUE))
  .stop_if_not(all(dhSegRows[-nrow(dhSegRows),2] < dhSegRows[-1,1], na.rm=TRUE))


  nbrOfSegments <- nrow(segs)
  verbose && cat(verbose, "Number of segments: ", nbrOfSegments)

  chromosome <- data$chromosome
  x <- data$x
  CT <- data$CT
  muN <- data$muN
  rho <- data$rho
  hasDH <- !is.null(rho)
  if (hasDH) {
    isHet <- !is.na(rho)
    isSnp <- isHet
  } else {
    isSnp <- !is.na(muN)
    isHet <- (isSnp & (muN == 1/2))
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Update the TCN segments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  chromosomes <- getChromosomes(fit)
  nbrOfChromosomes <- length(chromosomes)
  verbose && cat(verbose, "Number of chromosomes: ", nbrOfChromosomes)
  verbose && print(verbose, chromosomes)

  for (cc in seq_len(nbrOfChromosomes)) {
    chr <- chromosomes[cc]
    chrTag <- sprintf("chr%02d", chr)
    verbose && enter(verbose, sprintf("Chromosome %d ('%s') of %d", cc, chrTag, nbrOfChromosomes))
    rows <- which(is.element(segs[["chromosome"]], chr))
    verbose && cat(verbose, "Rows:")
    verbose && print(verbose, rows)

    segsCC <- segs[rows,,drop=FALSE]
    tcnSegRowsCC <- tcnSegRows[rows,,drop=FALSE]
    dhSegRowsCC <- dhSegRows[rows,,drop=FALSE]
    nbrOfSegmentsCC <- nrow(segsCC)
    verbose && cat(verbose, "Number of segments: ", nbrOfSegmentsCC)

    tcnIds <- sort(unique(segsCC[["tcnId"]]))
    I <- length(tcnIds)
    for (ii in seq_len(I)) {
      tcnId <- tcnIds[ii]
      verbose && enter(verbose, sprintf("TCN segment #%d ('%s') of %d", ii, tcnId, I))

      rowsII <- which(segsCC[["tcnId"]] == tcnId)
      J <- length(rowsII)
      # Nothing todo?
      if (!force && J == 1) {
        verbose && cat(verbose, "Nothing todo. Only one DH segmentation. Skipping.")
        verbose && exit(verbose)
        next
      }

      verbose && cat(verbose, "Rows:")
      verbose && print(verbose, rowsII)
      segsII <- segsCC[rowsII,,drop=FALSE]

      tcnSegRowsII <- tcnSegRowsCC[rowsII,,drop=FALSE]
      dhSegRowsII <- dhSegRowsCC[rowsII,,drop=FALSE]

      verbose && cat(verbose, "TCN & DH segRows before:")
      verbose && print(verbose, cbind(tcn=tcnSegRowsII, dh=dhSegRowsII))

      segRowsRange <- range(c(tcnSegRowsII, dhSegRowsII), na.rm=TRUE)
      verbose && printf(verbose, "Range [%d,%d]\n",
                                    segRowsRange[1], segRowsRange[2])

      tcnSegRowsIIBefore <- tcnSegRowsII
      nbrOfTCNsBefore <- segsII[1,"tcnNbrOfLoci"]
      # Sanity check
      .stop_if_not(diff(segRowsRange)+1L == nbrOfTCNsBefore)

      for (jj in seq_len(J)) {
        verbose && enter(verbose, sprintf("DH segment #%d of %d", jj, J))
        seg <- segsII[jj,,drop=FALSE]
        tcnSegRow <- unlist(tcnSegRowsII[jj,,drop=FALSE], use.names=FALSE)
        dhSegRow <- unlist(dhSegRowsII[jj,,drop=FALSE], use.names=FALSE)
        # Sanity check
        .stop_if_not(all(is.na(tcnSegRow)) || (tcnSegRow[1] <= tcnSegRow[2]))
        .stop_if_not(all(is.na(dhSegRow)) || (dhSegRow[1] <= dhSegRow[2]))

        # Sanity check
        idxsTCN <- tcnSegRow[1]:tcnSegRow[2]
        nbrOfTCNs <- sum(!is.na(CT[idxsTCN]))
        .stop_if_not(nbrOfTCNs == nbrOfTCNsBefore)

        if (joinSegments) {
          # (a) The TCN segment should have identical (start,end) boundaries as the DH region
          xStart <- seg[["dhStart"]]
          xEnd <- seg[["dhEnd"]]
          verbose && printf(verbose, "[xStart,xEnd] = [%.0f,%.0f]\n", xStart, xEnd)
          .stop_if_not(xStart <= xEnd)

          # (b) Identify units
          units <- which(chromosome == chr & xStart <= x & x <= xEnd)

          # (c) Drop units that are outside both the TCN and DH segments
          keep <- (segRowsRange[1] <= units & units <= segRowsRange[2])
          units <- units[keep]

          tcnSegRow <- range(units)
          verbose && printf(verbose, "[idxStart,idxEnd] = [%d,%d]\n", tcnSegRow[1], tcnSegRow[2])
          verbose && cat(verbose, "Number of TCN loci: ", length(units))

          # (c) Adjust for missing values
          keep <- which(!is.na(CT[units]))
          units <- units[keep]

          # (d) Adjust for DH boundaries
          if (jj > 1L) {
            minIdx <- tcnSegRowsII[jj-1L,2L, drop=TRUE]
            units <- units[units > minIdx]
          }
          if (jj < J) {
            maxIdx <- dhSegRowsII[jj+1L,1L, drop=TRUE]
            units <- units[units < maxIdx]
          }

          if (jj == J) {
#           maxIdx <- dhSegRowsII[jj+1L,1L, drop=TRUE]
#           units <- units[units < maxIdx]
          }

          tcnSegRow <- range(units)
          verbose && printf(verbose, "[idxStart,idxEnd] = [%d,%d]\n", tcnSegRow[1], tcnSegRow[2])
          verbose && cat(verbose, "Number of non-missing TCN loci: ", length(units))
        } else {
          throw("Not implemented yet.")  # /HB 2010-12-02
        } # if (joinSegments)

        gamma <- avgTCN(CT[units])
        # Sanity check
        .stop_if_not(length(units) == 0 || !is.na(gamma))

        # Update the segment boundaries, estimates and counts
        seg[["tcnStart"]] <- xStart
        seg[["tcnEnd"]] <- xEnd
        seg[["tcnMean"]] <- gamma
        seg[["tcnNbrOfLoci"]] <- length(units)
        seg[["tcnNbrOfSNPs"]] <- sum(isSnp[units])
        seg[["tcnNbrOfHets"]] <- sum(isHet[units])

        # Sanity check
        .stop_if_not(nrow(seg) == length(jj))

        segsII[jj,] <- seg
        tcnSegRowsII[jj,] <- tcnSegRow

        verbose && exit(verbose)
      } # for (jj ...)

      # Sanity check
      .stop_if_not(nrow(segsII) == length(rowsII))

      verbose && cat(verbose, "TCN & DH segRows afterward:")
      verbose && print(verbose, cbind(tcn=tcnSegRowsII, dh=dhSegRowsII))

##print(segsII)

      # Sanity check
      nbrOfTCNsAfter <- sum(segsII[,"tcnNbrOfLoci"], na.rm=TRUE)
      verbose && cat(verbose, "Number of TCNs before: ", nbrOfTCNsBefore)
      verbose && cat(verbose, "Number of TCNs after: ", nbrOfTCNsAfter)
      .stop_if_not(nbrOfTCNsAfter >= nbrOfTCNsBefore)

      # Sanity check
      .stop_if_not(nrow(dhSegRowsII) == nrow(tcnSegRowsII))
      .stop_if_not(all(tcnSegRowsII[,1] <= tcnSegRowsII[,2], na.rm=TRUE))
      .stop_if_not(all(tcnSegRowsII[-nrow(tcnSegRowsII),2] < tcnSegRowsII[-1,1], na.rm=TRUE))
      .stop_if_not(all(dhSegRowsII[,1] <= dhSegRowsII[,2], na.rm=TRUE))
      .stop_if_not(all(dhSegRowsII[-nrow(dhSegRowsII),2] < dhSegRowsII[-1,1], na.rm=TRUE))

      segsCC[rowsII,] <- segsII
      tcnSegRowsCC[rowsII,] <- tcnSegRowsII

      # Not needed anymore
      rowsII <- segsII <- NULL
      verbose && exit(verbose)
    } # for (ii ...)

    # Sanity check
    .stop_if_not(nrow(segsCC) == length(rows))

    # Sanity check
    .stop_if_not(nrow(dhSegRowsCC) == nrow(tcnSegRowsCC))
    .stop_if_not(all(tcnSegRowsCC[,1] <= tcnSegRowsCC[,2], na.rm=TRUE))
####################
if (!all(tcnSegRowsCC[-nrow(tcnSegRowsCC),2] < tcnSegRowsCC[-1,1], na.rm=TRUE)) {

  aa <- tcnSegRowsCC[-nrow(tcnSegRowsCC),2]
  bb <- tcnSegRowsCC[-1,1]
  delta <- bb - aa
  dd <- cbind(aa, bb, delta=delta)
  print(dd)
  dd <- subset(dd, delta == 0)
  print(dd)
  row <- dd[,1L,drop=TRUE]
  print(row)
  rr <- row + -10:10
  dd <- data[rr,]
  rownames(dd) <- rr
  print(dd)
print(tcnSegRowsII)
}
####################

    .stop_if_not(all(tcnSegRowsCC[-nrow(tcnSegRowsCC),2] < tcnSegRowsCC[-1,1], na.rm=TRUE))
    .stop_if_not(all(dhSegRowsCC[,1] <= dhSegRowsCC[,2], na.rm=TRUE))
    .stop_if_not(all(dhSegRowsCC[-nrow(dhSegRowsCC),2] < dhSegRowsCC[-1,1], na.rm=TRUE))

    segs[rows,] <- segsCC
    tcnSegRows[rows,] <- tcnSegRowsCC

    # Not needed anymore
    rows <- segsCC <- NULL
    verbose && exit(verbose)
  } # for (cc ...)

  # Sanity check
  .stop_if_not(nrow(dhSegRows) == nrow(tcnSegRows))
  .stop_if_not(all(tcnSegRows[,1] <= tcnSegRows[,2], na.rm=TRUE))
  .stop_if_not(all(tcnSegRows[-nrow(tcnSegRows),2] < tcnSegRows[-1,1], na.rm=TRUE))
  .stop_if_not(all(dhSegRows[,1] <= dhSegRows[,2], na.rm=TRUE))
  .stop_if_not(all(dhSegRows[-nrow(dhSegRows),2] < dhSegRows[-1,1], na.rm=TRUE))

  verbose && enter(verbose, "Update (C1,C2) per segment")
  # Append (C1,C2) estimates
  tcn <- segs$tcnMean
  dh <- segs$dhMean
  C1 <- 1/2*(1-dh)*tcn
  C2 <- tcn - C1
  segs$c1Mean <- C1
  segs$c2Mean <- C2
  verbose && exit(verbose)

  # Return results
  keep <- which(is.finite(fit$output$chromosome))
  fitS <- fit
  fitS$data <- data
  fitS$output[keep,] <- segs
  fitS$tcnSegRows[keep,] <- tcnSegRows

  # Sanity check
  tcnSegRows <- fitS$tcnSegRows
  dhSegRows <- fitS$dhSegRows
  .stop_if_not(nrow(dhSegRows) == nrow(tcnSegRows))
  .stop_if_not(all(tcnSegRows[,1] <= tcnSegRows[,2], na.rm=TRUE))
  .stop_if_not(all(tcnSegRows[-nrow(tcnSegRows),2] < tcnSegRows[-1,1], na.rm=TRUE))
  .stop_if_not(all(dhSegRows[,1] <= dhSegRows[,2], na.rm=TRUE))
  .stop_if_not(all(dhSegRows[-nrow(dhSegRows),2] < dhSegRows[-1,1], na.rm=TRUE))

  # Update 'flavor'
  fitS$params$flavor <- gsub(",", "&", flavor, fixed=TRUE)

  verbose && exit(verbose)

  fitS
}, protected=TRUE) # postsegmentTCN()
