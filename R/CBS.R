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
})


setMethodS3("as.data.frame", "CBS", function(x, ...) {
  getSegments(x, splitter=FALSE, ...);
})

setMethodS3("getSignalType", "CBS", function(fit, ...) {
  type <- fit$signalType;
  if (is.null(type)) type <- as.character(NA);
  type;
})

setMethodS3("signalType", "CBS", function(fit, ...) {
  getSignalType(fit);
})

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

  # Drop chromosome splitters?
  if (!splitters) {
    isSplitter <- lapply(segs[-1], FUN=is.na);
    isSplitter <- Reduce("&", isSplitter);
    segs <- segs[!isSplitter,];
  }

  if (nrow(segs) > 0) {
    segs$id <- getSampleName(fit);
  }

  segs;
}, private=TRUE)


###########################################################################/**
# @RdocMethod append
#
# @title "Appends one segmentation result to another"
#
# \description{
#   @get "title".
# }
# 
# @synopsis
#
# \arguments{
#  \item{x, other}{The two @see "CBS" objects to be combined.}
#  \item{other}{A @see "PSCBS" object.}
#  \item{addSplit}{If @TRUE, a "divider" is added between chromosomes.}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a @see "CBS" object of the same class as argument \code{x}.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################  
setMethodS3("append", "CBS", function(x, other, addSplit=TRUE, ...) {
  # To please R CMD check
  this <- x;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'other':
  other <- Arguments$getInstanceOf(other, class(this)[1]);
  for (field in c("data", "output")) {
    thisFF <- this[[field]];
    otherFF <- other[[field]];
    if (ncol(otherFF) != ncol(thisFF)) {
      throw(sprintf("Cannot merge %s objects. Arguments 'other' and 'this' has different number of columns in field '%s': %s != %s", class(this)[1], field, ncol(otherFF), ncol(thisFF)));
    }
  }

  # Argument 'addSplit':
  addSplit <- Arguments$getLogical(addSplit);


  # Allocate results
  res <- this;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Locus data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  data <- getLocusData(this);
  res$data <- rbind(data, getLocusData(other));


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Segmentation data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  indexOffset <- nrow(data);
  fields <- c("output", "segRows");
  for (field in fields[-1]) {
    other[[field]] <- other[[field]] + indexOffset;
  }

  splitter <- if (addSplit) NA else NULL;
  for (field in fields) {
    res[[field]] <- rbind(this[[field]], splitter, other[[field]]);
    rownames(res[[field]]) <- NULL;
  }

  # Sanity check
  ns <- sapply(res[fields], FUN=nrow);
  stopifnot(all(ns == ns[1]));

  res;
}) # append()



setMethodS3("writeLocusData", "CBS", function(fit, filename=sprintf("%s,byLocus.tsv", getSampleName(fit)), path=NULL, sep="\t", nbrOfDecimals=4L, addHeader=TRUE, createdBy=NULL, overwrite=FALSE, skip=FALSE, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Arguments 'filename' & 'path':
  pathname <- Arguments$getWritablePathname(filename, path=path, mustNotExist=(!overwrite && !skip));

  # Argument 'nbrOfDecimals':
  nbrOfDecimals <- Arguments$getInteger(nbrOfDecimals);


  # File already exists?
  if (isFile(pathname)) {
    # Skip?
    if (skip) {
      return(pathname);
    }

    # Overwrite!
    file.remove(pathname);
  }

  # Write to temporary file
  pathnameT <- pushTemporaryFile(pathname);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  data <- getLocusData(fit, ...);

  # Round of floating points
  if (!is.null(nbrOfDecimals)) {
    cols <- colnames(data);
    for (key in cols) {
      values <- data[[key]];
      if (is.double(values)) {
        values <- round(values, digits=nbrOfDecimals);
        data[[key]] <- values;
      }
    } # for (key ...)
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Build header
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (addHeader) {
    sigmaDelta <- estimateStandardDeviation(fit, method="diff");
#    sigmaResiduals <- estimateStandardDeviation(fit, method="res");

    createdOn <- format(Sys.time(), format="%Y-%m-%d %H:%M:%S %Z");
    hdr <- c(
      sampleName=getSampleName(fit),
      segmentationMethod=sprintf("segment() of %s", attr(fit, "pkgDetails")),
      nbrOfLoci=nbrOfLoci(fit),
      nbrOfSegments=nbrOfSegments(fit),
      joinSegments=fit$params$joinSegments,
      signalType=getSignalType(fit),
      sigmaDelta=sprintf("%.4f", sigmaDelta),
#      sigmaResiduals=sprintf("%.4f", sigmaResiduals),
      createdBy=createdBy,
      createdOn=createdOn,
      nbrOfDecimals=nbrOfDecimals,
      nbrOfColumns=ncol(data),
      columnNames=paste(colnames(data), collapse=", "),
      columnClasses=paste(sapply(data, FUN=function(x) class(x)[1]), collapse=", ")
    );
    bfr <- paste("# ", names(hdr), ": ", hdr, sep="");

    cat(file=pathnameT, bfr, sep="\n");
  } # if (addHeader)

  write.table(file=pathnameT, data, append=TRUE, quote=FALSE, sep=sep, 
                                          row.names=FALSE, col.names=TRUE);

  pathname <- popTemporaryFile(pathnameT);

  pathname;  
}) # writeLocusData()



setMethodS3("writeSegments", "CBS", function(fit, filename=sprintf("%s.tsv", getSampleName(fit)), path=NULL, addHeader=TRUE, createdBy=NULL, sep="\t", nbrOfDecimals=4L, splitters=FALSE, overwrite=FALSE, skip=FALSE, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Arguments 'filename' & 'path':
  pathname <- Arguments$getWritablePathname(filename, path=path, mustNotExist=(!overwrite && !skip));

  # Argument 'nbrOfDecimals':
  nbrOfDecimals <- Arguments$getInteger(nbrOfDecimals);


  # File already exists?
  if (isFile(pathname)) {
    # Skip?
    if (skip) {
      return(pathname);
    }

    # Overwrite!
    file.remove(pathname);
  }

  # Write to temporary file
  pathnameT <- pushTemporaryFile(pathname);


  sampleName <- getSampleName(fit);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  data <- getSegments(fit, splitters=splitters);

  # Round of floating points
  if (!is.null(nbrOfDecimals)) {
    cols <- c("start", "end");
    for (key in cols) {
      values <- data[[key]];
      if (is.double(values)) {
        values <- round(values, digits=0);
        data[[key]] <- values;
      }
    } # for (key ...)

    cols <- colnames(data);
    cols <- setdiff(cols, c("chromosome", "start", "end", "nbrOfLoci"));
    for (key in cols) {
      values <- data[[key]];
      if (is.double(values)) {
        values <- round(values, digits=nbrOfDecimals);
        data[[key]] <- values;
      }
    } # for (key ...) 
 }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Build header
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (addHeader) {
    sigmaDelta <- estimateStandardDeviation(fit, method="diff");
#    sigmaResiduals <- estimateStandardDeviation(fit, method="res");

    createdOn <- format(Sys.time(), format="%Y-%m-%d %H:%M:%S %Z");
    hdr <- c(
      sampleName=sampleName,
      segmentationMethod=sprintf("segment() of %s", attr(fit, "pkgDetails")),
      nbrOfLoci=nbrOfLoci(fit),
      nbrOfSegments=nbrOfSegments(fit),
      joinSegments=fit$params$joinSegments,
      signalType=getSignalType(fit),
      sigmaDelta=sprintf("%.4f", sigmaDelta),
#      sigmaResiduals=sprintf("%.4f", sigmaResiduals),
      createdBy=createdBy,
      createdOn=createdOn,
      nbrOfDecimals=nbrOfDecimals,
      nbrOfColumns=ncol(data),
      columnNames=paste(colnames(data), collapse=", "),
      columnClasses=paste(sapply(data, FUN=function(x) class(x)[1]), collapse=", ")
    );
    bfr <- paste("# ", names(hdr), ": ", hdr, sep="");

    cat(file=pathnameT, bfr, sep="\n");
  } # if (addHeader)

  write.table(file=pathnameT, data, append=TRUE, quote=FALSE, sep=sep, 
                                          row.names=FALSE, col.names=TRUE);

  pathname <- popTemporaryFile(pathnameT);

  pathname;  
}) # writeSegments()


############################################################################
# HISTORY:
# 2011-10-02
# o CLEANUP: Moved getChromosomes(), nbrOfChromosomes(), nbrOfSegments(),
#   nbrOfLoci() and print() to AbstractCBS.
# o Now the CBS class extends the AbstractCBS class.
# 2011-09-04
# o Added writeSegments() for CBS.
# o Added writeLocusData() for CBS.
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
