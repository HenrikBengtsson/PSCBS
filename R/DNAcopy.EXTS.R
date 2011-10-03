###########################################################################/**
# @set class=CBS
# @RdocMethod as.DNAcopy
#
# @title "Coerces a CBS object to a DNAcopy object"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{fit}{A @see "CBS" object."}
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns a @see "DNAcopy" object (of the \pkg{DNAcopy} package).
# }
#
# \examples{
#   @include "../incl/segmentByCBS.Rex"
#   @include "../incl/as.DNAcopy.Rex"
# }
#
# @author
#
# \seealso{
#   \code{\link[PSCBS:as.CBS.DNAcopy]{as.CBS()}}.
#   @seeclass
# }
#
# @keyword internal
#*/########################################################################### 
setMethodS3("as.DNAcopy", "CBS", function(fit, ...) {
  sampleName <- getSampleName(fit);
  if (is.na(sampleName)) sampleName <- "<NA>";

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup the 'data' field
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  data <- getLocusData(fit);

  # Keep only certain columns
  keep <- match(c("chromosome", "x"), colnames(data));
  keep <- c(keep, 3L);
  data <- data[,keep];

  # Sanity check
  stopifnot(ncol(data) == 3);

  # Rename column names
  colnames(data) <- c("chrom", "maploc", sampleName);

  class(data) <- c("CNA", "data.frame");
  attr(data, "data.type") <- "logratio";


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup the 'output' field
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  output <- getSegments(fit, splitter=FALSE);
  rownames <- rownames(output);

  output <- data.frame(
    ID        = sampleName,
    chrom     = output$chromosome,
    loc.start = output$start,
    loc.end   = output$end,
    num.mark  = output$nbrOfLoci,
    seg.mean  = output$mean,
    stringsAsFactors=FALSE
  );
  rownames(output) <- rownames;


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup up 'DNAcopy' object
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  res <- list();
  res$data <- data;
  res$output <- output;
  res$call <- NA;
  class(res) <- "DNAcopy";

  res;
}) # as.DNAcopy()



setMethodS3("nbrOfSegments", "DNAcopy", function(fit, ...) {
  segs <- fit$output;
  nrow(segs);
})

setMethodS3("nbrOfLoci", "DNAcopy", function(fit, ...) {
  nrow(fit$data);
})

setMethodS3("nbrOfSamples", "DNAcopy", function(fit, ...) {
  length(getSampleNames(fit, ...));
})

setMethodS3("getSampleNames", "DNAcopy", function(fit, ...) {
  names <- colnames(fit$data);
  names <- setdiff(names, c("chrom", "maploc"));
  names;
})

setMethodS3("getChromosomes", "DNAcopy", function(fit, ...) {
  chromosomes <- fit$data$chrom;
  sort(unique(chromosomes));
})


setMethodS3("estimateStandardDeviation", "DNAcopy", function(fit, sample=1, method=c("diff", "abs", "res"), estimator=c("mad", "sd"), na.rm=TRUE, weights=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'sample':
  sample <- Arguments$getIndex(sample, max=nbrOfSamples(fit));

  # Argument 'method':
  method <- match.arg(method);

  # Argument 'estimator':
  estimator <- match.arg(estimator);



  nbrOfLoci <- nbrOfLoci(fit);

  # Argument 'weights':
  if (!is.null(weights)) {
    weights <- Arguments$getNumerics(weights, range=c(0,Inf), length=rep(nbrOfLoci, times=2));
  }


  # Nothing to do?
  if (nbrOfLoci <= 1) {
    return(as.double(NA));
  }

  # Get the estimator function
  if (!is.null(weights)) {
    estimator <- sprintf("weighted %s", estimator);
    estimator <- R.utils::toCamelCase(estimator);
  }  
  estimatorFcn <- get(estimator, mode="function");


  # Extract sample of interest
  fit <- subset(fit, samplelist=sample);

  y <- fit$data[,3];

  if (method == "diff") {
    y <- diff(y);

    # Weighted estimator?
    if (!is.null(weights)) {
      # Calculate weights per pair
      weights <- (weights[1:(nbrOfLoci-1)]+weights[2:nbrOfLoci])/2;
      sigma <- estimatorFcn(y, w=weights, na.rm=na.rm)/sqrt(2);
    } else {
      sigma <- estimatorFcn(y, na.rm=na.rm)/sqrt(2);
    }
  } else if (method == "abs") {
    if (!is.null(weights)) {
      sigma <- estimatorFcn(y, w=weights, na.rm=na.rm);
    } else {
      sigma <- estimatorFcn(y, na.rm=na.rm);
    }
  } else if (method == "res") {
     yS <- extractSegmentMeansByLocus(fit);
     dy <- y - yS;
      if (!is.null(weights)) {
        sigma <- estimatorFcn(dy, w=weights, na.rm=na.rm);
      } else {
        sigma <- estimatorFcn(dy, na.rm=na.rm);
      }
  } else {
    throw("Method no implemented: ", method);
  }

  sigma;
}) # estimateStandardDeviation()


setMethodS3("extractSegmentMeansByLocus", "DNAcopy", function(fit, sample=1, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'sample':
  sample <- Arguments$getIndex(sample, max=nbrOfSamples(fit));

  # Extract sample of interest
  fit <- subset(fit, samplelist=sample);

  data <- fit$data;
  chr <- data$chrom;
  x <- data$maploc;
  y <- data[,3];

  segs <- fit$output;
  nbrOfSegments <- nrow(segs);
  nbrOfLoci <- nbrOfLoci(fit);

  yS <- y;
  for (ss in seq(length=nbrOfSegments)) {
    seg <- segs[ss,];
    idxs <- which(seg$chrom == chr & seg$loc.start <= x & x <= seg$loc.end);
    idxs <- Arguments$getIndices(idxs, max=nbrOfLoci);

    ySS <- y[idxs];
    ok <- is.finite(ySS);
    # Sanity check
    ## stopifnot(sum(ok) == seg$num.mark); # Not dealing with ties
    mu <- mean(ySS[ok]);
    yS[idxs] <- mu;
  } # for (ss ...)

  yS;
}, private=TRUE) # extractSegmentMeansByLocus()


############################################################################
# HISTORY:
# 2011-09-04
# o as.DNAcopy() did not drop "splitters" for the segment table.
# 2011-09-03
# o Added as.DNAcopy() for CBS to coerce a CBS object to a DNAcopy object.
# 2011-09-02
# o Added internal extractSegmentMeansByLocus() for DNAcopy, which is
#   used by estimateStandardDeviation(..., method="res").
# o Added estimateStandardDeviation() for DNAcopy.
# o ROBUSTNESS: Now getSampleNames() drops columns 'chrom' and 'maploc',
#   instead of assuming their positions.
# o ROBUSTNESS: Now nbrOfSamples() utilizes getSampleNames().
# o Added nbrOfSegments(), nbrOfLoci(), nbrOfSamples(), getSampleNames()
#   and getChromosomes() for DNAcopy.
# HISTORY FROM PRIVATE SCRIPTS:
# 2011-07-20
# o Added support for estimateStandardDeviation(..., method="res").
# o Added extractSegmentMeansByLocus().
# 2011-07-18
# o Added getSampleNames().
# o Added plotTracks() for DNAcopy.
# o Added nbrOfSegments(), nbrOfLoci() and nbrOfSamples().
# 2011-07-17
# o Added estimateStandardDeviation() to DNAcopy objects.
############################################################################
