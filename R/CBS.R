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
#  \item{...}{Not used.}
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
#   @see "segmentByCBS"
# }
#
# \author{Henrik Bengtsson}
#*/###########################################################################  
setConstructorS3("CBS", function(...) {
  extend(list(data=NULL, output=NULL), "CBS");
})


setMethodS3("print", "CBS", function(x, ...) {
  print(as.character(x));
})


setMethodS3("as.character", "CBS", function(x, ...) {
  # To please R CMD check
  fit <- x;

  s <- sprintf("%s:", class(fit)[1]);

  s <- c(s, sprintf("Sample name: %s", getSampleName(fit)));

  s <- c(s, sprintf("Number of segments: %d", nbrOfSegments(fit)));

  s <- c(s, sprintf("Number of loci: %d", nbrOfLoci(fit)));

  n <- getSegments(fit)$nbrOfLoci;
  q <- quantile(n, probs=c(0.00, 0.05, 0.25, 0.50, 0.75, 0.95, 1.00), na.rm=TRUE);
  qs <- sprintf("%g [%s]", q, names(q));
  s <- c(s, sprintf("Number of loci per segment: %s", paste(qs, collapse=", ")));

  chrs <- getChromosomes(fit);
  s <- c(s, sprintf("Chromosomes: [%d] %s", length(chrs), hpaste(chrs)));

  GenericSummary(s);
})


setMethodS3("as.data.frame", "CBS", function(x, ...) {
  getSegments(x, ...);
})

setMethodS3("getSampleName", "CBS", function(fit, ...) {
  name <- fit$sampleName;
  if (is.null(name)) {
    name <- as.character(NA);
  }
  name;
}, private=TRUE)

setMethodS3("sampleName", "CBS", function(fit, ...) {
  getSampleName(fit);
}, private=TRUE)

"sampleName<-" <- function(x, value) {
  UseMethod("sampleName<-");
}

setMethodS3("sampleName<-", "CBS", function(x, value) {
  fit <- x;

  # Argument 'value':
  value <- Arguments$getCharacter(value);

  fit$sampleName <- value;
  fit;
}, private=TRUE, addVarArgs=FALSE)


setMethodS3("getLocusData", "CBS", function(fit, ...) {
  data <- fit$data;
  data;
}, private=TRUE)

setMethodS3("getSegments", "CBS", function(fit, ...) {
  fit$output;
}, private=TRUE)

setMethodS3("nbrOfLoci", "CBS", function(fit, ...) {
  nrow(getLocusData(fit));
})

setMethodS3("nbrOfSegments", "CBS", function(fit, ...) {
  nrow(getSegments(fit));
})

setMethodS3("nbrOfChromosomes", "CBS", function(fit, ...) {
  length(getChromosomes(fit, ...));
})

setMethodS3("getChromosomes", "CBS", function(fit, ...) {
  data <- getLocusData(fit);
  chromosomes <- data$chromosome;
  chromosomes <- chromosomes[!is.na(chromosomes)];
  sort(unique(chromosomes));
})


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
  res$data <- rbind(this$data, other$data);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Segmentation data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  indexOffset <- nrow(this$data);
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



############################################################################
# HISTORY:
# 2011-09-04
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
