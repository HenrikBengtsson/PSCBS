###########################################################################/**
# @RdocClass PSCBS
#
# @title "The PSCBS class"
#
# \description{
#  @classhierarchy
#
#  A PSCBS is an object containing results from parent-specific copy-number
#  (PSCN) segmentation.
# }
#
# \usage{PSCBS(fit=list(), ...)}
#
# \arguments{
#   \item{fit}{A @list structure containing the PSCN segmentation results.}
#   \item{...}{Not used.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "HB"
#
# \seealso{
#   @see "PairedPSCBS".
# }
#*/###########################################################################
setConstructorS3("PSCBS", function(fit=list(), ...) {
  # Argument 'fit':
  if (!is.list(fit)) {
    throw("Argument 'fit' is not a list: ", class(fit)[1])
  }

  extend(AbstractCBS(fit, ...), "PSCBS")
})


setMethodS3("as.data.frame", "PSCBS", function(x, ...) {
  getSegments(x, splitters=TRUE, ...)
}, protected=TRUE)


setMethodS3("getLocusSignalNames", "PSCBS", function(fit, ...) {
  c("CT", "rho")
}, protected=TRUE)

setMethodS3("getSegmentTrackPrefixes", "PSCBS", function(fit, ...) {
  c("tcn", "dh")
}, protected=TRUE)


setMethodS3("getLocusData", "PSCBS", function(fit, indices=NULL, fields=c("asis"), ...) {
  # Argument 'indices':
  if (!is.null(indices)) {
    indices <- Arguments$getIndices(indices)
  }

  # Argument 'fields':
  fields <- match.arg(fields)

  data <- fit$data

  # Return requested indices
  if (!is.null(indices)) {
    # Map of final indices to current indices
    map <- match(indices, data$index)

    # Extract/expand...
    data <- data[map,]

    # Sanity check
    .stop_if_not(nrow(data) == length(indices))
  }

  data
}, protected=TRUE) # getLocusData()



setMethodS3("isSegmentSplitter", "PSCBS", function(fit, ...) {
  segs <- fit$output

  isSplitter <- lapply(segs[-1], FUN=is.na)
  isSplitter <- Reduce("&", isSplitter)

  isSplitter
}, protected=TRUE)


###########################################################################/**
# @RdocMethod getSegments
#
# @title "Gets the segments"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{simplify}{If @TRUE, redundant and intermediate information is dropped.}#  \item{splitters}{If @TRUE, "splitters" between chromosomes are
#     preserved, otherwise dropped.}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a SxK @data.frame, where S in the number of segments,
#   and K is the number of segment-specific fields.
# }
#
# @author "HB"
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getSegments", "PSCBS", function(fit, simplify=FALSE, splitters=TRUE, addGaps=FALSE, ...) {
  # Argument 'splitters':
  splitters <- Arguments$getLogical(splitters)

  segs <- fit$output

  # Drop chromosome splitters?
  if (!splitters) {
    isSplitter <- isSegmentSplitter(fit)
    segs <- segs[!isSplitter,]
  }

  # Add splitters for "gaps"...
  if (splitters && addGaps) {
    # Chromosome gaps
    n <- nrow(segs)
    chrs <- segs$chromosome
    gapsAfter <- which(diff(chrs) != 0L)
    gapsAfter <- gapsAfter[!is.na(chrs[gapsAfter])]
    nGaps <- length(gapsAfter)
    if (nGaps > 0L) {
      idxs <- seq_len(n)
      values <- rep(NA_integer_, times=nGaps)
      idxs <- insert(idxs, ats=gapsAfter+1L, values=values)
      segs <- segs[idxs,]
    }

    # Other gaps
    n <- nrow(segs)
    chrs <- segs$chromosome
    starts <- segs$tcnStart[-1L]
    ends <- segs$tcnEnd[-n]
    gapsAfter <- which(starts != ends)
    onSameChr <- (chrs[gapsAfter+1L] == chrs[gapsAfter] )
    gapsAfter <- gapsAfter[onSameChr]
    nGaps <- length(gapsAfter)
    if (nGaps > 0L) {
      idxs <- seq_len(n)
      values <- rep(NA_integer_, times=nGaps)
      idxs <- insert(idxs, ats=gapsAfter+1L, values=values)
      segs <- segs[idxs,]
    }
  }

##  if (nrow(segs) > 0) {
##    segs$id <- getSampleName(fit)
##  }

  if (simplify) {
    # If joinSegments was used (i.e. (start,end) are equal for TCN and DH)...
    if (fit$params$joinSegments) {
      # Sanity check
      .stop_if_not(all(segs$tcnStart == segs$dhStart, na.rm=TRUE))
      .stop_if_not(all(segs$tcnEnd == segs$dhEnd, na.rm=TRUE))

      names <- colnames(segs)
      keep <- !is.element(names, c("dhStart", "dhEnd"))
      segs <- segs[,keep]
      names <- colnames(segs)
      names[names == "tcnStart"] <- "start"
      names[names == "tcnEnd"] <- "end"
      colnames(segs) <- names
    }

    # Drop bootstrap columns, if any
    names <- colnames(segs)
    keep <- (regexpr("_[0-9]+(|[.][0-9]+)%$", names) == -1)
    segs <- segs[,keep]
  }

  segs
}, private=TRUE)



setMethodS3("getChangePoints", "PSCBS", function(fit, ...) {
  # Already available?
  cps <- fit$changepoints
  if (!is.null(cps)) return(cps)

  segs <- getSegments(fit, splitters=TRUE)
  tcn <- segs[["tcnMean"]]
  dh <- segs[["dhMean"]]
  C1 <- (1-dh) * tcn / 2
  C2 <- tcn - C1
  n <- length(tcn)

  # Calculate observed (alpha, radius, manhattan, dc1, dc2) data
  D1 <- C1[-n] - C1[-1L]
  D2 <- C2[-n] - C2[-1L]
  cps <- data.frame(
    alpha = atan2(D2, D1), # Changepoint angles in (0,2*pi)
    radius = sqrt(D2^2 + D1^2),
    manhattan = abs(D2) + abs(D1),
    d1 = D1,
    d2 = D2
  )

  cps
}, private=TRUE) # getChangePoints()


setMethodS3("normalizeTotalCNs", "PSCBS", function(fit, targetTCN=2, ...) {
  ## Fit using locus-level data
  data <- getLocusData(fit, ...)
  C <- data$CT
  .stop_if_not(!is.null(C))
  mu <- median(C, na.rm=TRUE)
  scale <- targetTCN / mu

  ## (a) Rescale locus-level data
  C <- scale * C
  data$CT <- C
  rm(list="C")
  fitN <- setLocusData(fit, data)

  ## (b) Rescale segment-level data
  segs <- getSegments(fit)
  fields <- colnames(segs)
  cnFields <- grep("^(tcn|c1|c2)", fields, value=TRUE)
  cnFields <- grep("(Id|Start|End|NbrOf)", cnFields, value=TRUE, invert=TRUE)
  for (field in cnFields) {
    segs[[field]] <- scale * segs[[field]]
  }
  fitN <- setSegments(fitN, segs)

  invisible(fitN)
})
