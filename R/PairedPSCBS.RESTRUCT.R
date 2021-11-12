setMethodS3("extractSegments", "PairedPSCBS", function(this, idxs, ..., verbose=FALSE) {
  fit <- this

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  updateSegRows <- function(segRows, idxs=NULL) {
    verbose && str(verbose, segRows)
    if (!is.null(idxs)) {
      segRows <- segRows[idxs,,drop=FALSE]
    }
#    verbose && cat(verbose, "Number of segments: ", nrow(segRows))
#    verbose && str(verbose, segRows)

    # Treat splitters separately
    isSplitter <- (is.na(segRows[,1]) & is.na(segRows[,2]))

    ns <- segRows[,2] - segRows[,1] + 1L
#    verbose && cat(verbose, "Number of loci per segment:")
#    verbose && str(verbose, ns)

    ns <- ns[!isSplitter]
    from <- c(1L, cumsum(ns)[-length(ns)]+1L)
    to <- from + (ns - 1L)
    segRows[!isSplitter,1] <- from
    segRows[!isSplitter,2] <- to
    verbose && str(verbose, segRows)

    # Sanity check
    ns2 <- segRows[,2] - segRows[,1] + 1L
    ns2 <- ns2[!isSplitter]
    .stop_if_not(all(ns2 == ns))

    segRows
  } # updateSegRows()


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'idxs':
  idxs <- Arguments$getIndices(idxs, max=nbrOfSegments(fit, splitters=TRUE))

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    pushState(verbose)
    on.exit(popState(verbose))
  }


  verbose && enter(verbose, "Extracting subset of segments")

  verbose && cat(verbose, "Number of segments: ", length(idxs))
  verbose && str(verbose, idxs)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract data and estimates
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  data <- getLocusData(fit)
  tcnSegRows <- fit$tcnSegRows
  dhSegRows <- fit$dhSegRows
  segs <- getSegments(fit)
  params <- fit$params

  # Sanity checks
  .stop_if_not(all(!is.na(data$chromosome) & !is.na(data$x)))
  .stop_if_not(length(tcnSegRows) == length(dhSegRows))

  # Sanity checks
  if (!params$joinSegments) {
    stop("Cannot extract subset of segments unless CNs are segmented using joinSegments=TRUE.")
  }

  if (params$flavor == "tcn,dh") {
    stop("NOT IMPLEMENTED: Extracting a subset of segments is not supported for flavor '", params$flavor, "'.")
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Subset segments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Update table of segments")
  segsT <- segs[idxs,,drop=FALSE]
  verbose && str(verbose, segsT)
  verbose && exit(verbose)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Subset data accordingly
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Update locus data")

  segRows <- tcnSegRows
  segRows <- segRows[idxs,,drop=FALSE]
  from <- segRows[[1]]
  to <- segRows[[2]]
  ok <- (!is.na(from) & !is.na(to))
  from <- from[ok]
  to <- to[ok]
  keep <- logical(nrow(data))
  for (rr in seq_along(from)) {
    keep[from[rr]:to[rr]] <- TRUE
  }
  keep <- which(keep)
  verbose && printf(verbose, "Identified %d (%.2f%%) of %d data rows:\n", length(keep), 100*length(keep)/nrow(data), nrow(data))
  verbose && str(verbose, keep)

  dataT <- data[keep,,drop=FALSE]
  verbose && str(verbose, dataT)

  verbose && exit(verbose)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Update 'segRows'
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Update 'segRows'")

  segRows <- updateSegRows(tcnSegRows, idxs=idxs)
  d <- tcnSegRows[idxs,] - segRows
  # Sanity check
  .stop_if_not(identical(d[,1], d[,2]))
  d <- d[,1]
  verbose && cat(verbose, "Row deltas:")
  verbose && str(verbose, d)

  tcnSegRows <- tcnSegRows[idxs,,drop=FALSE] - d
  verbose && str(verbose, tcnSegRows)
  # Sanity checks
  segRows <- tcnSegRows
  .stop_if_not(suppressWarnings(max(segRows, na.rm=TRUE)) <= nrow(dataT))
  drow <- segRows[-1,1] - segRows[-nrow(segRows),2]
  if (!all(is.na(drow) | (drow > 0))) {
    print(segRows)
    stop("INTERNAL ERROR: Generated 'tcnSegRows' is invalid, because it contains overlapping data chunks.")
  }

  dhSegRows <- dhSegRows[idxs,,drop=FALSE] - d
  verbose && str(verbose, dhSegRows)
  # Sanity checks
  segRows <- dhSegRows
  .stop_if_not(suppressWarnings(max(segRows, na.rm=TRUE)) <= nrow(dataT))
  drow <- segRows[-1,1] - segRows[-nrow(segRows),2]
  .stop_if_not(all(is.na(drow) | (drow > 0)))
  if (!all(is.na(drow) | (drow > 0))) {
    print(segRows)
    stop("INTERNAL ERROR: Generated 'dhSegRows' is invalid, because it contains overlapping data chunks.")
  }

  verbose && exit(verbose)


  # Create new object
  res <- fit
  res$data <- dataT
  res$output <- segsT
  res$tcnSegRows <- tcnSegRows
  res$dhSegRows <- dhSegRows

  verbose && exit(verbose)

  res
}, protected=TRUE) # extractSegments()



###########################################################################/**
# @set "class=PairedPSCBS"
# @RdocMethod mergeTwoSegments
#
# @title "Merge two neighboring segments"
#
# \description{
#   @get "title" by recalculating segment statistics.
# }
#
# @synopsis
#
# \arguments{
#  \item{left}{An @integer specifying the segments (left, left+1)
#    to be merged.}
#  \item{update}{If @TRUE, segment statistics are updated.}
#  \item{verbose}{A @logical or a @see "R.utils::Verbose" object.}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a @see "PairedPSCBS" with one less segment.
# }
#
# @author "HB"
#
# \seealso{
#   To drop regions (a connected set of segments) see \code{dropRegions()}.
#   @seeclass
# }
#*/###########################################################################
setMethodS3("mergeTwoSegments", "PairedPSCBS", function(this, left, update=TRUE, verbose=FALSE, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  nbrOfSegments <- nbrOfSegments(this, splitters=TRUE)
  # Argument 'left':
  left <- Arguments$getIndex(left, max=nbrOfSegments-1L)

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    pushState(verbose)
    on.exit(popState(verbose))
  }


  verbose && enter(verbose, "Merging two segments")
  verbose && printf(verbose, "Segments to be merged: %s & %s\n", left, left+1)
  verbose && cat(verbose, "Number of segments before merging: ", nbrOfSegments)
  verbose && cat(verbose, "Number of segments after merging: ", nbrOfSegments-1L)

  segs <- getSegments(this, splitters=TRUE)
  tcnSegRows <- this$tcnSegRows
  dhSegRows <- this$dhSegRows

  rows <- c(left,left+1)
  segsT <- segs[rows,,drop=FALSE]

  # Sanity check
  chrs <- segsT[["chromosome"]]
  if (chrs[1] != chrs[2]) {
    stop("Cannot merge segments that are on different chromosomes: ", chrs[1], " != ", chrs[2])
  }

  # Merge segments
  segT <- segsT[1,]
  fields <- colnames(segsT)

  # (chromosome, tcnId, dhId)
  idxsUsed <- 1:3

  # Starts
  idxs <- grep("Start$", fields)
  T <- as.matrix(segsT[,idxs,drop=FALSE])
  segT[,idxs] <- colMins(T, na.rm=TRUE, useNames=FALSE)
  idxsUsed <- c(idxsUsed, idxs)

  # Ends
  idxs <- grep("End$", fields)
  T <- as.matrix(segsT[,idxs,drop=FALSE])
  segT[,idxs] <- colMaxs(T, na.rm=TRUE, useNames=FALSE)
  idxsUsed <- c(idxsUsed, idxs)

  # Counts
  idxs <- grep("NbrOf", fields)
  segT[,idxs] <- colSums(segsT[,idxs,drop=FALSE])
  idxsUsed <- c(idxsUsed, idxs)

  # "Invalidate" remaining entries
  if (update) {
    idxsTodo <- setdiff(seq_along(fields), idxsUsed)
    segT[,idxsTodo] <- NA
  }

  # Update segment table
  segs[rows[1],] <- segT
  segs <- segs[-rows[2],]

  # Update 'segRows' tables
  segRows <- tcnSegRows
  segRows[rows[1],2] <- segRows[rows[2],2]
  segRows <- segRows[-rows[2],]
  tcnSegRows <- segRows

  segRows <- dhSegRows
  segRows[rows[1],2] <- segRows[rows[2],2]
  segRows <- segRows[-rows[2],]
  dhSegRows <- segRows

  # Create results object
  res <- this
  res$output <- segs
  res$tcnSegRows <- tcnSegRows
  res$dhSegRows <- dhSegRows

  # Update the segment statistics?
  if (update) {
    res <- updateMeans(res)
  }

  verbose && exit(verbose)

  res
}, private=TRUE)
