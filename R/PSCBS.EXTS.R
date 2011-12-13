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


setMethodS3("getLocusData", "PSCBS", function(fit, indices=NULL, ...) {
  # Argument 'indices':
  if (!is.null(indices)) {
    indices <- Arguments$getIndices(indices);
  }

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

  data;
}, protected=TRUE)


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
