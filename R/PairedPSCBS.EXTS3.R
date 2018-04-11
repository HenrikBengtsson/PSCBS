setMethodS3("extractLocusLevelC1C2", "PairedPSCBS", function(fit, ...) {
  # Extract locus-level data
  data <- getLocusData(fit)
  C <- data$CT
  rho <- data$rho

  # Swapped (C1,C2) <-> (C2,C1) for some segments?
  fields <- colnames(getSegments(fit))
  if (is.element("c1c2Swap", fields)) {
    # FIXME: When PSCBS is updated.
    # WORKAROUND: extractSegmentDataByLocus() in PSCBS v0.40.4 requires:
    # that fields "chromosome", "tcnStart", "tcnEnd" are always requested
    # /2014-03-21
    c1c2Swap <- extractSegmentDataByLocus(fit, fields=c("c1c2Swap",
                                     "chromosome", "tcnStart", "tcnEnd"))
    c1c2Swap <- c1c2Swap[["c1c2Swap"]]
    if (any(c1c2Swap)) {
      rho[c1c2Swap] <- -rho[c1c2Swap]
    }
  }

  C1 <- 1/2*(1-rho)*C
  C2 <- C - C1

  data.frame(C1=C1, C2=C2)
}, private=TRUE) # extractLocusLevelC1C2()


setMethodS3("extractLocusLevelTCN", "PairedPSCBS", function(fit, ...) {
  data <- getLocusData(fit)
  C <- data$CT
}, private=TRUE) # extractLocusLevelTCN()



setMethodS3("extractDhSegment", "PairedPSCBS", function(fit, idx, what=c("hets", "SNPs", "all"), ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'what':
  what <- match.arg(what)


  segs <- getSegments(fit, splitters=TRUE)
  .stop_if_not(!is.null(segs))
  nbrOfSegments <- nrow(segs)

  # Argument 'idx':
  idx <- Arguments$getIndex(idx, max=nbrOfSegments)

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    pushState(verbose)
    on.exit(popState(verbose))
  }



  verbose && enter(verbose, "Extracting a specific DH segment")

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract the data and segmentation results
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  data <- getLocusData(fit)
  .stop_if_not(!is.null(data))

  segs <- getSegments(fit, splitters=TRUE)
  .stop_if_not(!is.null(segs))

  verbose && enter(verbose, "Subsetting segment")
  # Subset the region-level data
  seg <- segs[idx,,drop=FALSE]

  isDivider <- all(is.na(seg))
  if (isDivider) {
    verbose && cat("Cannot extract DH segment. Not a valid segment: ", idx)
    verbose && exit(verbose)
    return(NULL)
  }

  verbose && print(verbose, seg)
  verbose && cat(verbose, "Number of TCN markers: ", sum(seg[["tcnNbrOfLoci"]], na.rm=TRUE))
  verbose && exit(verbose)

  verbose && enter(verbose, "Subsetting data")
  units <- seq_len(nrow(data))

  # Keep only chromosome of interest
  chr <- as.numeric(seg[,"chromosome"])
  if (!is.na(chr)) {
    keep <- which(data$chromosome == chr)
    units <- units[keep]
    data <- data[keep,]
  }

  # Keep only loci within the segment
  xRange <- as.numeric(seg[,c("dhStart", "dhEnd")])
  keep <- which(xRange[1] <= data$x & data$x <= xRange[2])
  units <- units[keep]
  data <- data[keep,]

  muN <- data$muN
  isSnp <- is.finite(muN)

  # Keep only SNPs?
  if (is.element(what, c("SNPs", "hets"))) {
    keep <- which(isSnp)
    units <- units[keep]
    data <- data[keep,]
  }

  # Keep only heterozygous SNPs?
  if (what == "hets") {
    isHet <- (muN == 1/2)
    keep <- which(isHet)
    units <- units[keep]
    data <- data[keep,]
  }
  verbose && exit(verbose)

  n <- nrow(data)
  verbose && cat(verbose, "Number of loci in DH segment: ", n)

  # Special case?
  listOfDhLociNotPartOfSegment <- fit$listOfDhLociNotPartOfSegment
  if (!is.null(listOfDhLociNotPartOfSegment)) {
    tcnId <- seg[,"tcnId"]
    dhId <- seg[,"dhId"]
    dhLociNotPartOfSegment <- listOfDhLociNotPartOfSegment[[tcnId]]
    if (!is.null(dhLociNotPartOfSegment)) {
      lociToExclude <- dhLociNotPartOfSegment[[dhId]]
      verbose && cat(verbose, "Excluding loci that belongs to a flanking segment: ", length(lociToExclude))
      drop <- match(lociToExclude, units)
      units <- units[-drop]
      data <- data[-drop,]
      n <- nrow(data)
    }
  }

  verbose && cat(verbose, "Number of units: ", n)
  verbose && cat(verbose, "Number of TCN markers: ", seg[,"tcnNbrOfLoci"])

  # Sanity check
  if (what == "hets" && n > 0) .stop_if_not(n == seg[,"dhNbrOfLoci"])

  fitS <- fit
  fitS$data <- data
  fitS$output <- seg

  verbose && exit(verbose)

  fitS
}, protected=TRUE) # extractDhSegment()
