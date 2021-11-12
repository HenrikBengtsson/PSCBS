setMethodS3("c", "PSCBS", function(..., addSplit = TRUE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  args <- list(...)

  ## Nothing todo?
  nargs <- length(args)
  if (nargs == 1) return(args[[1]])

  isNA <- function(x) is.logical(x) && length(x) == 1L && is.na(x)

  res <- args[[1]]
  fields <- c("output", "tcnSegRows", "dhSegRows")
  
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
        stop(sprintf("Cannot concatenate %s and %s objects, because they have different sets of columns in field %s: {%s} [n=%d] != {%s} [n=%d]", sQuote(class(res)[1]), sQuote(class(arg)[1]), sQuote(field), paste(sQuote(colnames(data)), collapse=", "), ncol(data), paste(sQuote(colnames(data_arg)), collapse=", "), ncol(data_arg)))
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


setMethodS3("extractChromosomes", "PSCBS", function(x, chromosomes, ...) {
  # To please R CMD check
  this <- x

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'chromosomes':
  disallow <- c("NaN", "Inf")
  chromosomes <- Arguments$getIntegers(chromosomes, range=c(0,Inf), disallow=disallow)
  .stop_if_not(all(is.element(chromosomes, getChromosomes(this))))

  # Always extract in order
  chromosomes <- unique(chromosomes)
  chromosomes <- sort(chromosomes)

  # Allocate results
  res <- this

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Locus data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  chromosome <- NULL; rm(list="chromosome") # To please R CMD check
  res$data <- subset(res$data, chromosome %in% chromosomes)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Segmentation data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify rows to subset
  rows <- which(is.element(res$output$chromosome, chromosomes))
  for (field in c("output", "tcnSegRows", "dhSegRows")) {
    res[[field]] <- res[[field]][rows,,drop=FALSE]
  }

  # Identify chromosome offsets
  chrStarts <- match(getChromosomes(this), this$data$chromosome)
  chrEnds <- c(chrStarts[-1]-1L, nrow(this$data))
  chrLengths <- chrEnds - chrStarts + 1L

  chrLengthsExcl <- chrLengths

  keep <- match(chromosomes, getChromosomes(this))
  chrLengthsExcl[keep] <- 0L
  cumChrLengthsExcl <- cumsum(chrLengthsExcl)

  shifts <- cumChrLengthsExcl[keep]
  .stop_if_not(all(is.finite(shifts)))

  # Adjust indices
  for (cc in seq_along(chromosomes)) {
    chromosome <- chromosomes[cc]
    shift <- shifts[cc]
    # Nothing to do?
    if (shift == 0) next
    for (field in c("tcnSegRows", "dhSegRows")) {
      segRows <- res[[field]]
      rows <- which(res$output$chromosome == chromosome)
      segRows[rows,] <- segRows[rows,] - shift
      res[[field]] <- segRows
    }
  }

  res
}, protected=TRUE)
