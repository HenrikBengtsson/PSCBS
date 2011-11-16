###########################################################################/**
# @RdocClass CBS
#
# @title "The CBS class"
#
# \description{
#   A CBS object holds results from the
#   Circular Binary Segmentation (CBS) method
#   for \emph{one} sample for one or more chromosomes.
#
#  @classhierarchy
# }
# 
# @synopsis
#
# \arguments{
#  \item{...}{Arguments passed to the constructor of @see "AbstractCBS".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \section{Difference to DNAcopy object}{
#   A CBS object is similar to DNAcopy objects with the major
#   difference that a CBS object holds only one sample, whereas
#   a DNAcopy object can hold more than one sample.
# }
#
# \section{See also}{
#  The @see "segmentByCBS" method returns an object of this class.
# }
#
# \author{Henrik Bengtsson}
#*/###########################################################################  
setConstructorS3("CBS", function(...) {
  extend(AbstractCBS(list(data=NULL, output=NULL), ...), "CBS");
})


setMethodS3("all.equal", "CBS", function(target, current, check.attributes=FALSE, ...) {
  # Compare class attributes
  res <- all.equal(class(target), class(current));
  if (!isTRUE(res)) {
    return(res);
  }

  # WORKAROUND: segmentByCBS() return getSegments(fit)$id without NA:s for
  # splitters, unless append() is used.
  # TO DO: Fix segmentByCBS() /HB 2011-10-08
  segs <- getSegments(target);

  isSplitter <- isSegmentSplitter(target);
  segs[isSplitter, "sampleName"] <- NA;
  target$output <- segs;

  segs <- getSegments(current);
  isSplitter <- isSegmentSplitter(current);
  segs[isSplitter, "sampleName"] <- NA;
  current$output <- segs;

  NextMethod("all.equal", target, current, check.attributes=check.attributes, ...);
}, protected=TRUE)



###########################################################################/**
# @RdocMethod as.data.frame
#
# @title "Gets the table of segments"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns a @data.frame, where each row corresponds to 
#   a unique segment.
# }
# 
# @author
#
# \seealso{
#   Utilizes @seemethod "getSegments".
#   @seeclass.
# }
#
# @keyword internal
#*/###########################################################################  
setMethodS3("as.character", "CBS", function(x, ...) {
  # To please R CMD check
  fit <- x;

  s <- sprintf("%s:", class(fit)[1]);

  s <- c(s, sprintf("Sample name: %s", getSampleName(fit)));

  s <- c(s, sprintf("Signal type: %s", getSignalType(fit)));

  s <- c(s, sprintf("Number of segments: %d", nbrOfSegments(fit)));

  s <- c(s, sprintf("Number of loci: %d", nbrOfLoci(fit)));

  n <- getSegments(fit)$nbrOfLoci;
  q <- quantile(n, probs=c(0.00, 0.05, 0.25, 0.50, 0.75, 0.95, 1.00), na.rm=TRUE);
  qs <- sprintf("%g [%s]", q, names(q));
  s <- c(s, sprintf("Number of loci per segment: %s", paste(qs, collapse=", ")));

  chrs <- getChromosomes(fit);
  s <- c(s, sprintf("Chromosomes: [%d] %s", length(chrs), hpaste(chrs)));

  s <- c(s, sprintf("Standard deviation: %g", estimateStandardDeviation(fit)));

  tt <- grep("Call$", colnames(getLocusData(fit)), value=TRUE);
  s <- c(s, sprintf("Locus calls: [%d] %s", length(tt), hpaste(tt)));

  segs <- getSegments(fit);
  callCols <- grep("Call$", colnames(segs), value=TRUE);
  callTypes <- gsub("Call$", "", callCols);
  s <- c(s, sprintf("Types of segment calls: [%d] %s", length(callTypes), hpaste(callTypes)));
  for (kk in seq(along=callCols)) {
    key <- callCols[kk];
    type <- callTypes[kk];
    n <- sum(segs[,key], na.rm=TRUE);
    if (type == "loss") {
      nC <- sum(isWholeChromosomeLost(fit));
    } else if (type == "gain") {
      nC <- sum(isWholeChromosomeGained(fit));
    } else {
      nC <- NA;
    }
    s <- c(s, sprintf("Number of chromosomes (segments) called '%s': %d (%d)", type, nC, n));
  }

  GenericSummary(s);
}, protected=TRUE)


setMethodS3("as.data.frame", "CBS", function(x, ...) {
  getSegments(x, splitter=FALSE, ...);
}, protected=TRUE)


setMethodS3("getSignalType", "CBS", function(fit, ...) {
  type <- fit$signalType;
  if (is.null(type)) type <- as.character(NA);
  type;
}, protected=TRUE)


setMethodS3("signalType", "CBS", function(fit, ...) {
  getSignalType(fit);
}, protected=TRUE)


"signalType<-" <- function(x, value) {
  UseMethod("signalType<-");
}

setMethodS3("signalType<-", "CBS", function(x, value) {
  fit <- x;

  # Argument 'value':
  value <- Arguments$getCharacter(value);

  fit$signalType <- value;
  fit;
}, private=TRUE, addVarArgs=FALSE)



setMethodS3("getLocusData", "CBS", function(fit, addCalls=NULL, ...) {
  # Argument 'addCalls':
  if (is.logical(addCalls)) {
    addCalls <- Arguments$getLogical(addCalls);
    if (!addCalls) {
      addCalls <- NULL;
    }
  } else {
    addCalls <- Arguments$getCharacters(addCalls);
  }

  data <- fit$data;

  # Append segment calls?
  if (length(addCalls) > 0) {
    callsL <- extractCallsByLocus(fit);
    if (is.character(addCalls)) {
      callsL <- callsL[,addCalls];
    }

    # Sanity check
    stopifnot(nrow(callsL) == nrow(data));

    data <- cbind(data, callsL);
  }

  data;
}, private=TRUE) # getLocusData()


setMethodS3("isSegmentSplitter", "CBS", function(fit, ...) {
  segs <- fit$output;

  isSplitter <- lapply(segs[-1], FUN=is.na);
  isSplitter <- Reduce("&", isSplitter);

  isSplitter;
}, protected=TRUE)


setMethodS3("getSegments", "CBS", function(fit, splitters=TRUE, ...) {
  # Argument 'splitters':
  splitters <- Arguments$getLogical(splitters);

  segs <- fit$output;

  isSplitter <- isSegmentSplitter(fit);

  # Add 'sampleName' column?
  if (nrow(segs) > 0) {
    sampleName <- rep(getSampleName(fit), times=nrow(segs));
    sampleName[isSplitter] <- as.character(NA);
    if (!is.element("sampleName", colnames(segs))) {
      segs <- cbind(sampleName=I(sampleName), segs);
    } else {
      segs[,"sampleName"] <- sampleName;
    }
  }

  # Drop chromosome splitters?
  if (!splitters) {
    segs <- segs[!isSplitter,];
  }

  segs;
}, private=TRUE)


setMethodS3("getSegmentSizes", "CBS", function(fit, by=c("length", "count"), ...) {
  by <- match.arg(by);

  data <- getSegments(fit, ...);
  if (by == "length") {
    res <- data[["end"]]-data[["start"]]+1L;
  } else if (by == "count") {
    res <- data[["nbrOfLoci"]];
  }
  res;
})


setMethodS3("updateBoundaries", "CBS", function(fit, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }
 
  verbose && enter(verbose, "Updating boundaries");
  verbose && cat(verbose, "Number of segments: ", 
                                  nbrOfSegments(fit, splitters=FALSE));

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Extract the data and segmentation results
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  data <- getLocusData(fit);
  segs <- getSegments(fit, splitters=TRUE);
  segRows <- fit$segRows;

  nbrOfSegments <- nrow(segs);
  chromosome <- data$chromosome;
  x <- data$x;
  y <- data$y;
  w <- data$w;
  hasWeights <- !is.null(w);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Update segments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  for (ss in seq(length=nbrOfSegments)) {
    verbose && enter(verbose, sprintf("Segment %d of %d", ss, nbrOfSegments));
    segRow <- segRows[ss,];
    seg <- segs[ss,];

    # A splitter - nothing todo?
    if (is.na(segRow[[1]]) && is.na(segRow[[2]])) {
      next;
    }

    # (a) Identify units (loci)
    units <- segRow[[1]]:segRow[[2]];
    verbose && cat(verbose, "Loci:");
    verbose && str(verbose, units);

    # (b) Extract signals
    ySS <- y[units];
    xSS <- x[units];
    cSS <- chromosome[units];
    if (hasWeights) {
      wSS <- w[units];
    }

    # (c) Drop missing values
    keep <- (!is.na(ySS) & !is.na(xSS) & !is.na(cSS));
    if (hasWeights) {
      keep <- keep & (!is.na(wSS) & wSS > 0);
    }
    keep <- which(keep);
    ySS <- ySS[keep];
    xSS <- xSS[keep];
    cSS <- cSS[keep];
    if (hasWeights) {
      wSS <- wSS[keep];
    }
    units <- units[keep];
    verbose && cat(verbose, "Loci (non-missing):");
    verbose && str(verbose, units);

    # (d) Identify (chromosome, start, stop)
    stopifnot(all(cSS == cSS[1]));
    cSS <- cSS[1];
    xRange <- range(xSS, na.rm=TRUE);
    verbose && cat(verbose, "Range:");
    verbose && print(verbose, xRange);

    # (e) Update segment information
    seg$chromosome <- cSS;
    seg$start <- xRange[1];
    seg$end <- xRange[2];

    segs[ss,] <- seg;

    verbose && exit(verbose);
  } # for (ss ...)

  # Update results
  res <- fit;
  res$output <- segs;

  # Rejoin segments?
  if (isTRUE(res$params$joinSegments)) {
    res <- joinSegments(res, verbose=less(verbose,10));
  }

  verbose && exit(verbose);

  res;
}) # updateBoundaries()



setMethodS3("updateMeans", "CBS", function(fit, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }
 
  verbose && enter(verbose, "Updating mean level estimates");
  verbose && cat(verbose, "Number of segments: ", 
                                  nbrOfSegments(fit, splitters=FALSE));

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Extract the data and segmentation results
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  data <- getLocusData(fit);
  segs <- getSegments(fit, splitters=TRUE);
  segRows <- fit$segRows;

  nbrOfSegments <- nrow(segs);
  chromosome <- data$chromosome;
  x <- data$x;
  y <- data$y;
  w <- data$w;
  hasWeights <- !is.null(w);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Update segments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  for (ss in seq(length=nbrOfSegments)) {
    verbose && enter(verbose, sprintf("Segment %d of %d", ss, nbrOfSegments));
    segRow <- segRows[ss,];
    seg <- segs[ss,];

    # A splitter - nothing todo?
    if (is.na(segRow[[1]]) && is.na(segRow[[2]])) {
      next;
    }

    # (a) Identify units (loci)
    units <- segRow[[1]]:segRow[[2]];

    # (b) Extract signals
    ySS <- y[units];
    if (hasWeights) {
      wSS <- w[units];
    }

    # (c) Drop missing values
    keep <- (!is.na(ySS));
    if (hasWeights) {
      keep <- keep & (!is.na(wSS) & wSS > 0);
    }
    keep <- which(keep);
    ySS <- ySS[keep];
    if (hasWeights) {
      wSS <- wSS[keep];
    }
    units <- units[keep];
    nbrOfLoci <- length(units);

    # (d) Update mean
    if (hasWeights) {
      wSS <- wSS / sum(wSS);
      gamma <- sum(wSS*ySS);
    } else {
      gamma <- mean(ySS);
    }

    # Sanity check
    stopifnot(nbrOfLoci == 0 || !is.na(gamma));

    # (d) Update the segment statistics
    seg$mean <- gamma;
    seg$nbrOfLoci <- nbrOfLoci;

    segs[ss,] <- seg;

    verbose && exit(verbose);
  } # for (ss ...)

  # Return results
  res <- fit;
  res$output <- segs;

  verbose && exit(verbose);

  res;
}, protected=TRUE) # updateMeans()



############################################################################
# HISTORY:
# 2011-11-15
# o Now updateMeans() uses locus-specific weights, iff available.
# o Added updateBoundaries() for CBS to update (start,stop) per segment.
# o CORRECTNESS: Now updateMeans() for CBS identify loci by the internal
#   'segRows' field and no longer by locations of segment boundaries,
#   which gave slightly incorrect estimates for "tied" loci.
# 2011-10-16
# o Added isSegmentSplitter().
# 2011-10-08
# o Relabelled column 'id' to 'sampleName' returned by getSegments().
# o BUG FIX: getSegments() for CBS would not set 'id' for "splitter" rows.
# o Added mergeTwoSegments() for CBS.
# o Added updateMeans() for CBS.
# o Added all.equal() for CBS.
# 2011-10-02
# o CLEANUP: Moved getChromosomes(), nbrOfChromosomes(), nbrOfSegments(),
#   nbrOfLoci() and print() to AbstractCBS.
# o Now the CBS class extends the AbstractCBS class.
# 2011-09-04
# o Added getSignalType() for CBS.
# o Added argument 'addCalls' to getLocusData().
# o Added getSampleName() for CBS.
# 2011-09-03
# o Added print() and as.character() for CBS.
# o Added CBS() constructor. Although it rairly will be used
#   it we be a place holder for the documentation.
# 2011-09-02
# o Added nbrOfLoci(), nbrOfSegments(), nbrOfChromosomes() and
#   getChromosomes() for CBS.
# 2010-11-19
# o Added append() for CBS objects.
############################################################################
