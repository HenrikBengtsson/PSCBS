###########################################################################/**
# @set "class=CBS"
# @RdocMethod callGainsAndLosses
#
# @title "Calls gains and losses" 
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{adjust}{A positive scale factor adjusting the sensitivity of the
#    caller, where a value less (greater) than 1.0 makes the caller 
#    less (more) sensitive.}
#  \item{method}{A @character string specifying the calling algorithm to use.}
#  \item{...}{Additional/optional arguments used to override the default
#    parameters used by the caller.}
# }
#
# \value{
#  Returns a @see "PSCBS::CBS" object where @logical columns
#  'lossCall' and 'gainCall' have been appended to the segmentation table.
# }
#
# \section{The UCSF caller}{
#   "Previously, we defined a segment to be gained/lost if it were at 
#    least two times the sample MAD away from the median segmented value."
#    (D. Albertson 2011-07-18)
# }
#
# \examples{
#   @include "../incl/segmentByCBS.Rex"
#   @include "../incl/segmentByCBS,calls.Rex"
# }
# 
# @author
#
# \seealso{
#   @seemethod "callAmplifications".
#   @seemethod "callOutliers".
#   @seeclass
# }
#
# @keyword internal 
#*/###########################################################################  
setMethodS3("callGainsAndLosses", "CBS", function(fit, adjust=1.0, method=c("ucsf-mad"), ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'adjust':
  adjust <- Arguments$getDouble(adjust, range=c(0, Inf));

  # Argument 'method':
  method <- match.arg(method);

  userArgs <- list(...);

  params <- list();

  # Allocate calls
  naValue <- as.logical(NA);
  nbrOfSegments <- nbrOfSegments(fit);
  gainCalls <- lossCalls <- rep(naValue, times=nbrOfSegments);

  if (method == "ucsf-mad") {
    # Default arguments
    args <- list(
      chromosomes = intersect(getChromosomes(fit), c(0L, 1:22)),
      scale = 2.0
    );

    # Override by (optional) user-specified arguments
    for (key in names(userArgs)) {
      args[[key]] <- userArgs[[key]];
    }

    # Extract arguments
    chromosomes <- args$chromosomes;
    scale <- args$scale;

    # Argument check
    chromosomes <- Arguments$getVector(chromosomes, lengths=c(1,Inf));
    scale <- Arguments$getDouble(scale, range=c(0,Inf));

    # Estimate the whole-genome standard deviation of the TCNs
    sigma <- estimateStandardDeviation(fit, chromosomes=chromosomes, 
                                       method="res", estimator="mad");

    # Sanity check
    sigma <- Arguments$getDouble(sigma, range=c(0,Inf));

    # Calculate the threshold
    tau <- scale * sigma;

    # Make more or less sensitive
    tau <- adjust * tau;

    # Calculate segment statistics (utilizing DNAcopy methods)
    stats <- segments.summary(as.DNAcopy(fit));
    kept <- is.element(rownames(fit$output), rownames(stats));

    naValue <- as.double(NA);
    mu <- rep(naValue, times=nbrOfSegments);

    # The segmented mean levels
    mu[kept] <- stats$seg.median;

    # The median segmented level
    muR <- median(mu, na.rm=TRUE);

    # The threshold for losses
    tauLoss <- muR - tau;

    # The threshold for gains
    tauGain <- muR + tau;

    # Call
    lossCalls <- (mu <= tauLoss);   # Losses
    gainCalls <- (mu >= tauGain);   # Gains

    # Call parameters used
    params$method <- method;
    params$adjust <- adjust;
    params$sigmaMAD <- sigma;
    params$scale <- scale;
    params$muR <- muR;
    params$tau <- tau;
    params$tauLoss <- tauLoss;
    params$tauGain <- tauGain;
  }

  # Sanity check
  stopifnot(length(lossCalls) == nbrOfSegments);
  stopifnot(length(gainCalls) == nbrOfSegments);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Update 'DNAcopy' object
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # (a) segmentation table
  segs <- fit$output;
  segs$lossCall <- lossCalls;
  segs$gainCall <- gainCalls;
  fit$output <- segs;

  # (b) parameters
  allParams <- fit$params;
  if (is.null(allParams)) {
    allParams <- list();
  }
  allParams$callGainsAndLosses <- params;
  fit$params <- allParams;


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Return the updated 'CBS' object.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  fit;
}, private=TRUE) # callGainsAndLosses()





###########################################################################/**
# @set "class=CBS"
# @RdocMethod callAmplifications
#
# @title "Calls (focal) amplifications" 
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{adjust}{A positive scale factor adjusting the sensitivity of the
#    caller, where a value less (greater) than 1.0 makes the caller 
#    less (more) sensitive.}
#  \item{maxLength}{A @double scalar specifying the maximum length of a segment
#    in order for it to be considered a focal amplification.}
#  \item{method}{A @character string specifying the calling algorithm to use.}
#  \item{...}{Additional/optional arguments used to override the default
#    parameters used by the caller.}
# }
#
# \value{
#  Returns a @see "PSCBS::CBS" object where @logical column 
#  'amplificationCall' has been appended to the segmentation table.
# }
#
# \section{The UCSF caller}{
#   "Previously, we defined a segment to be gained/lost if it were at 
#    least two times the sample MAD away from the median segmented value."
#    (D. Albertson 2011-07-18)
# }
#
# @author
#
# \references{
#   [1] Fridlyand et al. (2006)\cr
# }
#
# \seealso{
#   @seemethod "callGainsAndLosses".
#   @seemethod "callOutliers".
#   @seeclass
# }
#
# @keyword internal 
#*/###########################################################################  
setMethodS3("callAmplifications", "CBS", function(fit, adjust=1.0, maxLength=20e6, method=c("ucsf-exp"), ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'adjust':
  adjust <- Arguments$getDouble(adjust, range=c(0, Inf));

  # Argument 'maxLength':
  maxLength <- Arguments$getDouble(maxLength, range=c(0, Inf));

  # Argument 'method':
  method <- match.arg(method);


  userArgs <- list(...);

  params <- list();

  # Allocate calls
  naValue <- as.logical(NA);
  nbrOfSegments <- nbrOfSegments(fit);
  calls <- rep(naValue, times=nbrOfSegments);

  if (method == "ucsf-exp") {
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # Call arguments
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # Default arguments
    args <- list(
      minLevel = 0.0,
      lambda   = 1.0,
      degree   = 3
    );

    # Override by (optional) user-specified arguments
    for (key in names(userArgs)) {
      args[[key]] <- userArgs[[key]];
    }

    # Extract arguments
    minLevel <- args$minLevel;
    lambda <- args$lambda;
    degree <- args$degree;

    # Validate arguments
    minLevel <- Arguments$getDouble(minLevel, range=c(-Inf, Inf));
    lambda <- Arguments$getDouble(lambda, range=c(0, Inf));
    degree <- Arguments$getDouble(degree, range=c(1, Inf));


    segs <- fit$output;

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # Rule #1: Only consider segments that are short enough
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # The lengths (in bp) of the segments
    start <- segs$start;
    end <- segs$end;
    length <- end - start + 1L;
    keep1 <- (length <= maxLength);

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # Rule #2: Only consider segments that have a mean level
    #          that is large enough.
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # The mean levels of the segments
    mu <- segs$mean;
    keep2 <- (mu >= minLevel);

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # Rule #3: Only consider segments that have a mean level
    #          that is much larger than either of the 
    #          flanking segments.
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # The mean levels of the flanking segments
    muL <- c(NA, mu[-nbrOfSegments]);
    muR <- c(mu[-1], NA);

    # The difference in mean levels to the flanking segments
    deltaL <- mu - muL;
    deltaR <- mu - muR;

    # The maximum difference to either of the flanking segments
    delta <- pmax(deltaL, deltaR, na.rm=TRUE);

    # The threshold for calling segments amplified
    tau <- exp(-lambda * mu^degree);

    # Make more or less sensitive
    tau <- adjust * tau;

    keep3 <- (delta >= tau);

    # Amplification calls
    calls <- (keep1 & keep2 & keep3);

    # Call parameters used
    params$method <- method;
    params$adjust <- adjust;
    params$maxLength <- maxLength;
    params$minLevel <- minLevel;
    params$lambda <- lambda;
    params$degree <- degree;
    params$tau <- tau;
  }

  # Sanity check
  stopifnot(length(calls) == nbrOfSegments);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Update 'DNAcopy' object
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # (a) segmentation table
  segs <- fit$output;
  segs$amplificationCall <- calls;
  fit$output <- segs;

  # (b) parameters
  allParams <- fit$params;
  if (is.null(allParams)) {
    allParams <- list();
  }
  allParams$callAmplifications <- params;
  fit$params <- allParams;


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Return the updated 'CBS' object.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  fit;
}, private=TRUE) # callAmplifications()



###########################################################################/**
# @set "class=CBS"
# @RdocMethod callOutliers
#
# @title "Calls outliers" 
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{adjust}{A positive scale factor adjusting the sensitivity of the
#    caller, where a value less (greater) than 1.0 makes the caller 
#    less (more) sensitive.}
#  \item{method}{A @character string specifying the calling algorithm to use.}
#  \item{...}{Additional/optional arguments used to override the default
#    parameters used by the caller.}
# }
#
# \value{
#  Returns a @see "PSCBS::CBS" object where @logical column
#  'outlierCall' has been appended to the segmentation table.
# }
#
# \section{The UCSF caller}{
#  "Finally, to identify single technical or biological outliers such 
#   as high level amplifications, the presence of the outliers within 
#   a segment was allowed by assigning the original observed log2ratio
#   to the clones for which the observed values were more than four
#   tumor-specific MAD away from the smoothed values." [1; Suppl. Mat.]
# }
#
# @author
#
# \references{
#   [1] Fridlyand et al. (2006)\cr
# }
#
# \seealso{
#   @seemethod "callGainsAndLosses".
#   @seemethod "callAmplifications".
#   @seeclass
# }
#
# @keyword internal 
#*/###########################################################################  
setMethodS3("callOutliers", "CBS", function(fit, adjust=1.0, method=c("ucsf-mad"), ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'adjust':
  adjust <- Arguments$getDouble(adjust, range=c(0, Inf));

  # Argument 'method':
  method <- match.arg(method);


  userArgs <- list(...);

  params <- list();

  # Allocate calls
  nbrOfLoci <- nbrOfLoci(fit);
  naValue <- as.logical(NA);
  calls <- rep(naValue, times=nbrOfLoci);

  if (method == "ucsf-mad") {
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # Call arguments
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # Default arguments
    args <- list(
      scale = 4.0
    );

    # Override by (optional) user-specified arguments
    for (key in names(userArgs)) {
      args[[key]] <- userArgs[[key]];
    }

    # Extract arguments
    scale <- args$scale;

    # Validate arguments
    scale <- Arguments$getDouble(scale, range=c(0, Inf));


    # Genomic annotations
    data <- fit$data;
    chromosome <- data$chromosome;
    x <- data$x;

    # CN signals
    data <- fit$data;
    y <- data[,3];

    # Segmented CN signals
    yS <- extractSegmentMeansByLocus(fit);

    # CN residuals (relative to segment means)
    dy <- y - yS;

    segs <- fit$output;

    # Allocate per-segment SD estimates
    nbrOfSegments <- nbrOfSegments(fit);
    naValue <- as.double(NA);
    sds <- rep(naValue, times=nbrOfSegments);

    naValue <- as.double(NA);
    for (ss in seq(length=nbrOfSegments)) {
      seg <- segs[ss,];

      # Identify loci in current segment
      idxs <- which(seg$chromosome == chromosome & 
                    seg$start <= x & x <= seg$end);

      # Sanity check
      idxs <- Arguments$getIndices(idxs, max=nbrOfLoci);

      # Extract CN residuals
      dySS <- dy[idxs];

      # Calculate MAD for segment
      sdSS <- mad(dySS, na.rm=TRUE);

      # Threshold for outliers
      tau <- scale * sdSS;

      # Make more or less sensitive
      tau <- adjust * tau;

      # Call outliers
      naValue <- as.logical(NA);
      callsSS <- rep(naValue, times=length(dySS));
      callsSS[-tau <= dySS & dySS <= +tau] <- FALSE;
      callsSS[dySS > +tau] <- TRUE;
      callsSS[dySS < -tau] <- TRUE;

      # Record
      calls[idxs] <- callsSS;
      sds[ss] <- sdSS;
    } # for (ss ...)

    params$method <- method;
    params$adjust <- adjust;
    params$scale <- scale;
    params$sds <- sds;
  }


  # Sanity check
  stopifnot(length(calls) == nbrOfLoci);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Update 'DNAcopy' object
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # (a) segmentation table
  data <- fit$data;
  data$outlierCall <- calls;
  fit$data <- data;

  # (b) parameters
  allParams <- fit$params;
  if (is.null(allParams)) {
    allParams <- list();
  }
  allParams$callOutliers <- params;
  fit$params <- allParams;


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Return the updated 'CBS' object.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  fit;
}, private=TRUE) # callOutliers()





############################################################################
# HISTORY:
# 2011-09-04
# o Adopted the calling methods from ditto of the DNAcopy class.
# 2011-09-01
# o Now callGainsAndLosses() returns a DNAcopy where the segmentation
#   table has the new column 'tcnCall'.
# 2011-08-19
# o Added argument 'callParams' to plotTracks() for DNAcopy.
# 2011-07-24
# o Added callOutliers().
# 2011-07-21
# o Now amplified segments are also highlighted.
# 2011-07-20
# o Added callAmplifications().
# 2011-07-20
# o Now callGainsAndLosses() estimates the noise level on autosomes only.
# o Now callGainsAndLosses() returns parameters used.
# o Updated callGainsAndLosses() to estimate the std. dev. as the
#   MAD of the *residuals* (not the absolute) values.
# o Added support for estimateStandardDeviation(..., method="res").
# o Added extractSegmentMeansByLocus().
# o Added drawCentromeres().
# 2011-07-18
# o Added getSampleNames().
# o Added plotTracks() for DNAcopy.
# o Added callGainsAndLosses() to DNAcopy objects.
# o Added nbrOfSegments(), nbrOfLoci() and nbrOfSamples().
# 2011-07-17
# o Added estimateStandardDeviation() to DNAcopy objects.
############################################################################
