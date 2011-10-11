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
  isSplitter <- lapply(segs[-1], FUN=is.na);
  isSplitter <- Reduce("&", isSplitter);
  segs[isSplitter, "sampleName"] <- NA;
  target$output <- segs;

  segs <- getSegments(current);
  isSplitter <- lapply(segs[-1], FUN=is.na);
  isSplitter <- Reduce("&", isSplitter);
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


setMethodS3("getSegments", "CBS", function(fit, splitters=TRUE, ...) {
  # Argument 'splitters':
  splitters <- Arguments$getLogical(splitters);

  segs <- fit$output;

  isSplitter <- lapply(segs[-1], FUN=is.na);
  isSplitter <- Reduce("&", isSplitter);

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

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Extract the data and segmentation results
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  data <- getLocusData(fit);

  segs <- getSegments(fit);
  keep <- is.finite(segs$chromosome);
  segs <- segs[keep,,drop=FALSE];

  nbrOfSegments <- nrow(segs);
  verbose && cat(verbose, "Number of segments: ", nbrOfSegments);

  chromosome <- data$chromosome;
  x <- data$x;
  y <- data$y;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Update segments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  for (ss in seq(length=nbrOfSegments)) {
    verbose && enter(verbose, sprintf("Segment %d of %d", ss, nbrOfSegments));
    seg <- segs[ss,];

    chr <- seg[["chromosome"]];
    chrTag <- sprintf("chr%02d", chr);

    keys <- c("start", "end");
    xStart <- seg[[keys[1]]];
    xEnd <- seg[[keys[2]]];
    verbose && printf(verbose, "[xStart,xEnd] = [%.0f,%.0f]\n", xStart, xEnd);
    # Nothing todo?
    if (is.na(xStart) && is.na(xEnd)) {
      next;
    }

    stopifnot(xStart <= xEnd);

    # (b) Identify units
    units <- which(chromosome == chr & xStart <= x & x <= xEnd);

    # (c) Adjust for missing values
    value <- y;
    keep <- which(!is.na(value[units]));
    units <- units[keep];

    # (d) Update mean
    gamma <- mean(value[units]);

    # Sanity check
    stopifnot(length(units) == 0 || !is.na(gamma));

    # Update the segment boundaries, estimates and counts
    key <- "mean";
    seg[[key]] <- gamma;

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
