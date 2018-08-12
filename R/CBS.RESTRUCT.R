setMethodS3("shiftTCN", "CBS", function(fit, shift, update=TRUE, ...) {
  # Argument 'shift':
  shift <- Arguments$getDouble(shift, disallow=c("NA", "NaN", "Inf"))

  data <- getLocusData(fit)
  data$y <- data$y + shift
  fit$data <- data
  # Not needed anymore
  data <- NULL

  if (update) {
    fit <- updateMeans(fit, ...)
  }

  fit
}, protected=TRUE)


###########################################################################/**
# @set "class=CBS"
# @RdocMethod c
# @alias c.PSCBS
#
# @title "Concatenates segmentation results"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{\dots}{One or more @see "AbstractCBS" objects to be combined.}
#  \item{addSplit}{If @TRUE, a "divider" is added between chromosomes.}
# }
#
# \value{
#   Returns an @see "AbstractCBS" object of the same class in \dots.
# }
#
# @author "HB"
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("c", "CBS", function(..., addSplit = TRUE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  args <- list(...)

  ## Nothing todo?
  nargs <- length(args)
  if (nargs == 1) return(args[[1]])

  isNA <- function(x) is.logical(x) && length(x) == 1L && is.na(x)

  res <- args[[1]]
  fields <- c("output", "segRows")
  
  for (ii in 2:nargs) {
    arg <- args[[ii]]

    if (isNA(arg)) {
      if (addSplit) {
        warning(sprintf("Detected explicit NA in call to c(<%s>, ..., addSplit = TRUE). Ignoring", class(args[[1]])[1]))
        next
      }
      ## Add "splitter"
      for (field in fields) {
        res[[field]] <- rbind(res[[field]], NA)
      }
    } else {
      ## Locus-level data
      data <- getLocusData(res)
      data_arg <- getLocusData(arg)
      if (!all(colnames(data_arg) == colnames(data))) {
        throw(sprintf("Cannot concatenate %s and %s objects, because they have different sets of columns in field %s: {%s} [n=%d] != {%s} [n=%d]", sQuote(class(res)[1]), sQuote(class(arg)[1]), sQuote(field), paste(sQuote(colnames(data)), collapse=", "), ncol(data), paste(sQuote(colnames(data_arg)), collapse=", "), ncol(data_arg)))
      }

      indexOffset <- nrow(data)
      
      data <- rbind(data, getLocusData(arg))
      res[["data"]] <- data
      
      # Segmentation data
      for (field in fields[-1]) {
        arg[[field]] <- arg[[field]] + indexOffset
      }
      splitter <- if (addSplit) NA else NULL
      for (field in fields) {
        res[[field]] <- rbind(res[[field]], splitter, arg[[field]])
      }

      # Known segments
      ksT <- res$params$knownSegments
      ksT$length <- NULL  # In case it's been added
      ksO <- arg$params$knownSegments
      ksO$length <- NULL  # In case it's been added
      res$params$knownSegments <- rbind(ksT, ksO)
    }
  } ## for (ii ...)

  ## Drop row names, iff they've been added
  for (field in fields) rownames(res[[field]]) <- NULL
  
  # Sanity check
  ns <- sapply(res[fields], FUN = nrow)
  .stop_if_not(all(ns == ns[1]))

  res
}) # c()



setMethodS3("extractSegments", "CBS", function(this, idxs, ..., verbose=FALSE) {
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
  segRows <- fit$segRows
  segs <- getSegments(fit)
  params <- fit$params

  # Sanity checks
  .stop_if_not(all(!is.na(data$chromosome) & !is.na(data$x)))

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

  segRowsT <- segRows[idxs,,drop=FALSE]
  from <- segRowsT[[1]]
  to <- segRowsT[[2]]
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
  segRowsT <- updateSegRows(segRowsT)
  d <- segRows[idxs,] - segRowsT

  # Sanity check
  .stop_if_not(identical(d[,1], d[,2]))
  d <- d[,1]
  verbose && cat(verbose, "Row deltas:")
  verbose && str(verbose, d)

  segRows <- segRows[idxs,,drop=FALSE] - d
  verbose && str(verbose, segRows)
  # Sanity checks
  .stop_if_not(suppressWarnings(max(segRows, na.rm=TRUE)) <= nrow(dataT))
  drow <- segRows[-1,1] - segRows[-nrow(segRows),2]
  .stop_if_not(all(is.na(drow) | (drow > 0)))
  if (!all(is.na(drow) | (drow > 0))) {
    print(segRows)
    throw("INTERNAL ERROR: Generated 'segRows' is invalid, because it contains overlapping data chunks.")
  }

  verbose && exit(verbose)


  # Create new object
  res <- fit
  res$data <- dataT
  res$output <- segsT
  res$segRows <- segRows

  verbose && exit(verbose)

  res
}, protected=TRUE) # extractSegments()



setMethodS3("mergeTwoSegments", "CBS", function(this, left, update=TRUE, verbose=FALSE, ...) {
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

  segs <- getSegments(this)
  segRows <- this$segRows

  rows <- c(left,left+1)
  segsT <- segs[rows,,drop=FALSE]

  # Sanity check
  chrs <- segsT[["chromosome"]]
  if (chrs[1] != chrs[2]) {
    throw("Cannot merge segments that are on different chromosomes: ", chrs[1], " != ", chrs[2])
  }

  # Merge segments
  segT <- segsT[1,]
  fields <- colnames(segsT)
  idxsUsed <- c()

  # (id) [as in label]
  idxs <- grep("(I|i)d$", fields)
  idxsUsed <- c(idxsUsed, idxs)

  # (chromosome)
  idxs <- grep("chromosome$", fields)
  idxsUsed <- c(idxsUsed, idxs)

  # Starts
  idxs <- grep("(S|s)tart$", fields)
  T <- as.matrix(segsT[,idxs,drop=FALSE])
  segT[,idxs] <- colMins(T, na.rm=TRUE)
  idxsUsed <- c(idxsUsed, idxs)

  # Ends
  idxs <- grep("(E|e)nd$", fields)
  T <- as.matrix(segsT[,idxs,drop=FALSE])
  segT[,idxs] <- colMaxs(T, na.rm=TRUE)
  idxsUsed <- c(idxsUsed, idxs)

  # Counts
  idxs <- grep("(N|n)brOf", fields)
  segT[,idxs] <- colSums(segsT[,idxs,drop=FALSE])
  idxsUsed <- c(idxsUsed, idxs)

  # "Invalidate" remaining entries
  idxsTodo <- setdiff(seq_along(fields), idxsUsed)
  segT[,idxsTodo] <- NA

  # Update segment table
  segs[rows[1],] <- segT
  segs <- segs[-rows[2],]

  # Update 'segRows' tables
  segRows[rows[1],2] <- segRows[rows[2],2]
  segRows <- segRows[-rows[2],]

  # Create results object
  res <- this
  res$output <- segs
  res$segRows <- segRows

  # Update the segment statistics?
  if (update) {
    res <- updateMeans(res)
  }

  verbose && exit(verbose)

  res
}, protected=TRUE) # mergeTwoSegments()
