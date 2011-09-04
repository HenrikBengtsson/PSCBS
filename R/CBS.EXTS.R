###########################################################################/**
# @set class=DNAcopy
# @RdocMethod as.CBS
#
# @title "Coerces a DNAcopy object to a CBS object"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{fit}{A @see "DNAcopy" object (of the \pkg{DNAcopy} package.)}
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns a @see "CBS" object.
# }
#
# @author
#
# \seealso{
#   \code{\link[PSCBS:as.DNAcopy.CBS]{as.DNAcopy()}}.
#   @seeclass
# }
#
# @keyword internal
#*/########################################################################### 
setMethodS3("as.CBS", "DNAcopy", function(fit, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup the 'data' field
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  data <- fit$data;

  # Rename column names
  names <- colnames(data);
  names[names == "chrom"] <- "chromosome";
  names[names == "maploc"] <- "x";
  colnames(data) <- names;

  # Drop 'CNA' class and DNAcopy attributes
  class(data) <- c("data.frame");
  attr(data, "data.type") <- NULL;


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup the 'output' field
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  output <- fit$output;

  # Rename column names
  names <- colnames(output);
  names[names == "ID"] <- "id";
  names[names == "chrom"] <- "chromosome";
  names[names == "loc.start"] <- "start";
  names[names == "loc.end"] <- "end";
  names[names == "num.mark"] <- "nbrOfLoci";
  names[names == "seg.mean"] <- "mean";
  colnames(output) <- names;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup up 'CBS' object
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  res <- list();
  res$data <- data;
  res$output <- output;
  res$params <- list();
  class(res) <- c("CBS");

  res;
}) # as.CBS()



setMethodS3("extractByChromosomes", "CBS", function(x, chromosomes, ...) {
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
  data <- this$data;
  class <- class(data);
  class(data) <- "data.frame";
  data <- subset(data, chromosome %in% chromosomes);
  class(data) <- class;
  res$data <- data;


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Segmentation data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify rows to subset
  rows <- which(is.element(res$output$chromosome, chromosomes));
  for (field in c("output", "segRows")) {
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
    for (field in c("segRows")) {
      segRows <- res[[field]];
      rows <- which(res$output$chromosome == chromosome);
      segRows[rows,] <- segRows[rows,] - shift;
      res[[field]] <- segRows;
    }
  }

  res;
})


setMethodS3("subset", "CBS", function(x, chromlist=NULL, ...) {
  extractByChromosomes(x, chromosomes=chromlist, ...);
}, private=TRUE)



###########################################################################/**
# @set "class=CBS"
# @RdocMethod extractSegmentMeansByLocus
#
# @title "Extracts segments means at each locus" 
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#  Returns a @numeric @vector of length @seemethod "nbrOfLoci".
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @keyword internal 
#*/###########################################################################  
setMethodS3("extractSegmentMeansByLocus", "CBS", function(fit, ...) {
  data <- fit$data;
  chromosome <- data$chromosome;
  x <- data$x;
  y <- data[,3];

  segs <- fit$output;
  nbrOfSegments <- nrow(segs);
  nbrOfLoci <- nbrOfLoci(fit);

  yS <- y;
  for (ss in seq(length=nbrOfSegments)) {
    seg <- segs[ss,];
    idxs <- which(seg$chromosome == chromosome & 
                  seg$start <= x & x <= seg$end);
    idxs <- Arguments$getIndices(idxs, max=nbrOfLoci);

    ySS <- y[idxs];
    ok <- is.finite(ySS);

    # Sanity check
    ## stopifnot(sum(ok) == seg$nbrOfLoci); # Not dealing with ties

    mu <- mean(ySS[ok]);
    yS[idxs] <- mu;
  } # for (ss ...)

  yS;
}, private=TRUE) # extractSegmentMeansByLocus()



###########################################################################/**
# @set "class=CBS"
# @RdocMethod estimateStandardDeviation
#
# @title "Estimates the whole-genome standard deviation of the signals"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{chromosomes}{An optional @vector specifying the subset of 
#    chromosomes used for the estimate.  If @NULL, all chromosomes are used.}
#  \item{method}{A @character string specifying the method used.}
#  \item{estimator}{A @character string or a @function specifying the
#    internal estimator.}
#  \item{na.rm}{If @TRUE, missing values are dropped, otherwise not.}
#  \item{weights}{An optional @double @vector of @seemethod "nbrOfLoci" 
#    non-negative weights.}
#  \item{...}{Not used.}
# }
#
# \value{
#  Returns a non-negative @numeric scale.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @keyword internal 
#*/###########################################################################  
setMethodS3("estimateStandardDeviation", "CBS", function(fit, chromosomes=NULL, method=c("diff", "res", "abs"), estimator=c("mad", "sd"), na.rm=TRUE, weights=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'chromosomes':
  if (!is.null(chromosomes)) {
  }
 
  # Argument 'method':
  method <- match.arg(method);

  # Argument 'estimator':
  estimator <- match.arg(estimator);


  # Argument 'weights':
  if (!is.null(weights)) {
    nbrOfLoci <- nbrOfLoci(fit);
    weights <- Arguments$getNumerics(weights, range=c(0,Inf), 
                                     length=rep(nbrOfLoci, times=2));
  }


  # Get the estimator function
  if (!is.null(weights)) {
    estimator <- sprintf("weighted %s", estimator);
    estimator <- R.utils::toCamelCase(estimator);
  }  
  estimatorFcn <- get(estimator, mode="function");


  # Subset by chromosomes?
  if (!is.null(chromosomes)) {
    fit <- extractByChromosomes(fit, chromosomes=chromosomes);
  }

  nbrOfLoci <- nbrOfLoci(fit);
  # Nothing to do?
  if (nbrOfLoci <= 1) {
    sigma <- as.double(NA);
    attr(sigma, "nbrOfLoci") <- nbrOfLoci;
    attr(sigma, "df") <- as.integer(NA);
    return(sigma);
  }

  y <- fit$data[,3];

  if (method == "diff") {
    dy <- diff(y);
    # Weighted estimator?
    if (!is.null(weights)) {
      # Calculate weights per pair
      weights <- (weights[1:(nbrOfLoci-1)]+weights[2:nbrOfLoci])/2;
      sigma <- estimatorFcn(dy, w=weights, na.rm=na.rm)/sqrt(2);
    } else {
      sigma <- estimatorFcn(dy, na.rm=na.rm)/sqrt(2);
    }
    df <- length(dy);
  } else if (method == "res") {
    yS <- extractSegmentMeansByLocus(fit);
    dy <- y - yS;
    if (!is.null(weights)) {
      sigma <- estimatorFcn(dy, w=weights, na.rm=na.rm);
    } else {
      sigma <- estimatorFcn(dy, na.rm=na.rm);
    }
    df <- length(dy);
  } else if (method == "abs") {
    if (!is.null(weights)) {
      sigma <- estimatorFcn(y, w=weights, na.rm=na.rm);
    } else {
      sigma <- estimatorFcn(y, na.rm=na.rm);
    }
    df <- length(y);
  } else {
    throw("Method no implemented: ", method);
  }

  attr(sigma, "nbrOfLoci") <- nbrOfLoci;
  attr(sigma, "df") <- df;

  sigma;
}) # estimateStandardDeviation()



############################################################################
# HISTORY:
# 2011-09-04
# o Added estimateStandardDeviation() for CBS.
# o Added extractSegmentMeansByLocus() for CBS.
# 2011-09-03
# o Added as.CBS() for DNAcopy to coerce a DNAcopy object to a CBS object.
# 2011-09-02
# o Added extractByChromosomes() for CBS.
# o Added subset() for CBS for backward compatibility.
# o Added nbrOfLoci(), getChromosomes() and getSampleNames() for CBS.
# 2010-11-19
# o Added append() for CBS objects.
############################################################################
