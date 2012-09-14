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
# @author
#
# \seealso{
#   @see "PairedPSCBS".
# }
#*/###########################################################################
setConstructorS3("PSCBS", function(fit=list(), ...) {
  # Argument 'fit':
  if (!is.list(fit)) {
    throw("Argument 'fit' is not a list: ", class(fit)[1]);
  }

  extend(AbstractCBS(fit, ...), "PSCBS");
})


setMethodS3("as.data.frame", "PSCBS", function(x, ...) {
  getSegments(x, splitter=TRUE, ...);
}, protected=TRUE)


setMethodS3("getLocusSignalNames", "PSCBS", function(fit, ...) {
  c("CT", "rho");
}, protected=TRUE)

setMethodS3("getSegmentTrackPrefixes", "PSCBS", function(fit, ...) {
  c("tcn", "dh");
}, protected=TRUE)


setMethodS3("getLocusData", "PSCBS", function(fit, indices=NULL, fields=c("asis", "full"), ...) {
  # Argument 'indices':
  if (!is.null(indices)) {
    indices <- Arguments$getIndices(indices);
  }

  # Argument 'fields':
  fields <- match.arg(fields);


  data <- fit$data;

  # Return requested indices
  if (!is.null(indices)) {
    # Map of final indices to current indices
    map <- match(indices, data$index);

    # Extract/expand...
    data <- data[map,];

    # Sanity check
    stopifnot(nrow(data) == length(indices));
  }

  if (fields == "full") {
    names <- colnames(data);

    # Genotype calls
    if (!is.element("muN", names)) {
      require("aroma.light") || throw("Package not loaded: aroma.light");
      data$muN <- callNaiveGenotypes(data$betaN);
    }
    data$isHet <- (data$muN == 1/2);

    data$rho <- 2*abs(data$betaT-1/2);
    data$c1 <- 1/2*(1-data$rho)*data$CT;
    data$c2 <- data$CT - data$c1;

    # TumorBoost BAFs
    if (!is.element("betaTN", names)) {
      require("aroma.light") || throw("Package not loaded: aroma.light");
      data$betaTN <- normalizeTumorBoost(betaN=data$betaN, betaT=data$betaT, muN=data$muN);
    }
    data$rhoN <- 2*abs(data$betaTN-1/2);
    data$c1N <- 1/2*(1-data$rhoN)*data$CT;
    data$c2N <- data$CT - data$c1N;

    data$isSNP <- (!is.na(data$betaT) | !is.na(data$betaN));
    data$type <- ifelse(data$isSNP, "SNP", "non-polymorphic locus");

    # Labels
    data$muNx <- c("AA", "AB", "BB")[2*data$muN + 1L];
    data$isHetx <- c("AA|BB", "AB")[data$isHet + 1L];
  }

  data;
}, protected=TRUE) # getLocusData()


setMethodS3("isSegmentSplitter", "PSCBS", function(fit, ...) {
  segs <- fit$output;

  isSplitter <- lapply(segs[-1], FUN=is.na);
  isSplitter <- Reduce("&", isSplitter);

  isSplitter;
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
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################  
setMethodS3("getSegments", "PSCBS", function(fit, simplify=FALSE, splitters=TRUE, ...) {
  # Argument 'splitters':
  splitters <- Arguments$getLogical(splitters);

  segs <- fit$output;

  # Drop chromosome splitters?
  if (!splitters) {
    isSplitter <- isSegmentSplitter(fit);
    segs <- segs[!isSplitter,];
  }

##  if (nrow(segs) > 0) {
##    segs$id <- getSampleName(fit);
##  }

  if (simplify) {
    # If joinSegments was used (i.e. (start,end) are equal for TCN and DH)...
    if (fit$params$joinSegments) {
      # Sanity check
      stopifnot(all(segs$tcnStart == segs$dhStart, na.rm=TRUE));
      stopifnot(all(segs$tcnEnd == segs$dhEnd, na.rm=TRUE));

      names <- colnames(segs);
      keep <- !is.element(names, c("dhStart", "dhEnd"));
      segs <- segs[,keep];
      names <- colnames(segs);
      names[names == "tcnStart"] <- "start";
      names[names == "tcnEnd"] <- "end";
      colnames(segs) <- names;
    }

    # Drop bootstrap columns, if any
    names <- colnames(segs);
    keep <- (regexpr("_[0-9]+(|[.][0-9]+)%$", names) == -1);
    segs <- segs[,keep];
  }

  segs;
}, private=TRUE) 



############################################################################
# HISTORY:
# 2012-04-21
# o CLEANUP: Moved getSegmentSizes() from PSCBS to AbstractCBS.
# 2012-04-21
# o CLEANUP: Moved getSegmentSizes() from PairedPSCBS to PSCBS.
# 2012-02-27
# o Added argument 'fields' to getLocusData() for PairedPSCBS.
# 2011-12-12
# o Added optional argument 'indices' to getLocusData() to be able
#   to retrieve the locus-level data as indexed by input data.
# 2011-12-03
# o Added argument 'simplify' to getSegments().
# 2011-10-16
# o Added isSegmentSplitter().
# 2011-10-02
# o Now the CBS class extends the AbstractCBS class.
# o Added print() and as.data.frame() to PSCBS.
# o Added getSegments() to PSCBS.
# o DOCUMENTATION: Added Rdoc for several PSCBS methods.
# o Added a PSCBS constructor with documentation.
# 2010-12-01
# o Now also extractByChromosomes() and append() for PSCBS recognizes
#   fields 'tcnLociToExclude' and 'dhLociToExclude'.
# o BUG FIX: extractByChromosome() for PSCBS would call it self instead
#   of extractByChromosomes().
# 2010-11-26
# o Added extractByChromosomes() for PSCBS.
# 2010-09-26
# o getChromosomes() no longer returns NA divers.
# 2010-09-24
# o Added append() and more for PSCBS objects.
############################################################################
