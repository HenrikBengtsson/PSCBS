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
})


setMethodS3("getLocusData", "PSCBS", function(fit, ...) {
  data <- fit$data;
  data;
})

 
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
#  \item{splitters}{If @TRUE, "splitters" between chromosomes are 
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
setMethodS3("getSegments", "PSCBS", function(fit, splitters=TRUE, ...) {
  # Argument 'splitters':
  splitters <- Arguments$getLogical(splitters);

  segs <- fit$output;

  # Drop chromosome splitters?
  if (!splitters) {
    isSplitter <- lapply(segs[-1], FUN=is.na);
    isSplitter <- Reduce("&", isSplitter);
    segs <- segs[!isSplitter,];
  }

##  if (nrow(segs) > 0) {
##    segs$id <- getSampleName(fit);
##  }

  segs;
}, private=TRUE) 



setMethodS3("extractByChromosomes", "PSCBS", function(x, chromosomes, ...) {
  # To please R CMD check
  this <- x;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'chromosomes':
  chromosomes <- Arguments$getIntegers(chromosomes, disallow=c("NaN", "Inf"));
  stopifnot(all(is.element(chromosomes, getChromosomes(this))));

  # Always extract in order
  chromosomes <- unique(chromosomes);
  chromosomes <- sort(chromosomes);

  # Allocate results
  res <- this;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Locus data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  chromosome <- NULL; rm(chromosome); # To please R CMD check
  res$data <- subset(res$data, chromosome %in% chromosomes);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Segmentation data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify rows to subset
  rows <- which(is.element(res$output$chromosome, chromosomes));
  for (field in c("output", "tcnSegRows", "dhSegRows")) {
    res[[field]] <- res[[field]][rows,,drop=FALSE];
  }

  # Identify chromosome offsets
  chrStarts <- match(getChromosomes(this), this$data$chromosome);
  chrEnds <- c(chrStarts[-1]-1L, nrow(this$data));
  chrLengths <- chrEnds - chrStarts + 1L;

  chrLengthsExcl <- chrLengths;

  keep <- match(chromosomes, getChromosomes(this));
  chrLengthsExcl[keep] <- 0L;
  cumChrLengthsExcl <- cumsum(chrLengthsExcl);

  shifts <- cumChrLengthsExcl[keep];
  stopifnot(all(is.finite(shifts)));

  # Adjust indices
  for (cc in seq(along=chromosomes)) {
    chromosome <- chromosomes[cc];
    shift <- shifts[cc];
    # Nothing to do?
    if (shift == 0) next;
    for (field in c("tcnSegRows", "dhSegRows")) {
      segRows <- res[[field]];
      rows <- which(res$output$chromosome == chromosome);
      segRows[rows,] <- segRows[rows,] - shift;
      res[[field]] <- segRows;
    }
  }

  res;
})



###########################################################################/**
# @set "class=PSCBS"
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
#  \item{x, other}{The two @see "PSCBS" objects to be combined.}
#  \item{other}{A @see "PSCBS" object.}
#  \item{addSplit}{If @TRUE, a "divider" is added between chromosomes.}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a @see "PSCBS" object of the same class as argument \code{x}.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################  
setMethodS3("append", "PSCBS", function(x, other, addSplit=TRUE, ...) {
  # To please R CMD check
  this <- x;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'other':
  other <- Arguments$getInstanceOf(other, "PSCBS");
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
  res$data <- rbind(this$data, other$data);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Segmentation data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  indexOffset <- nrow(this$data);
  fields <- c("output", "tcnSegRows", "dhSegRows");
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



############################################################################
# HISTORY:
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
