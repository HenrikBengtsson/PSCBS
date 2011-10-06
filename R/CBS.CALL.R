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
#   If \code{method == "ucsf-mad"}, then segments are called using [1], i.e.
#   a segment is called gained or lost if its segment level is
#   at least two standard deviations away from the median segment level
#   on Chr1-22, where standard deviation is estimated using MAD.
# }
#
# \examples{
#   @include "../incl/segmentByCBS.Rex"
#   @include "../incl/segmentByCBS,calls.Rex"
# }
# 
# @author
#
# \references{
#   [1] Fridlyand et al. \emph{Breast tumor copy number aberration 
#       phenotypes and genomic instability}, BMC Cancer, 2006. \cr
# }
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
  segs <- getSegments(fit, splitter=TRUE);
  nbrOfRows <- nrow(segs);
  gainCalls <- lossCalls <- rep(naValue, times=nbrOfRows);

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
    kept <- is.element(rownames(segs), rownames(stats));

    naValue <- as.double(NA);
    mu <- rep(naValue, times=nbrOfRows);

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
  stopifnot(length(lossCalls) == nbrOfRows);
  stopifnot(length(gainCalls) == nbrOfRows);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Update 'DNAcopy' object
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # (a) segmentation table
  segs <- getSegments(fit, splitter=TRUE);
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
#   If \code{method == "ucsf-exp"}, then segments are called using [1], i.e.
#   a segment is called an amplification if ...
# }
#
# @author
#
# \references{
#   [1] Fridlyand et al. \emph{Breast tumor copy number aberration 
#       phenotypes and genomic instability}, BMC Cancer, 2006. \cr
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
#  Returns a @see "PSCBS::CBS" object where @logical columns
#  'negOutlierCall' and 'posOutlierCall' have been appended
#  to the segmentation table.
# }
#
# \section{The UCSF caller}{
#   If \code{method == "ucsf-mad"}, then loci are called using [1];
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
#   [1] Fridlyand et al. \emph{Breast tumor copy number aberration 
#       phenotypes and genomic instability}, BMC Cancer, 2006. \cr
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
  negOutlierCall <- posOutlierCall <- rep(naValue, times=nbrOfLoci);

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
    data <- getLocusData(fit);
    chromosome <- data$chromosome;
    x <- data$x;

    # CN signals
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
      callsSS[-tau <= dySS & dySS <= +tau] <- 0L;
      callsSS[dySS > +tau] <- +1L;
      callsSS[dySS < -tau] <- -1L;

      # Record
      negOutlierCall[idxs] <- (callsSS < 0L);
      posOutlierCall[idxs] <- (callsSS > 0L);

      sds[ss] <- sdSS;
    } # for (ss ...)

    params$method <- method;
    params$adjust <- adjust;
    params$scale <- scale;
    params$sds <- sds;
  }


  # Sanity check
  stopifnot(length(negOutlierCall) == nbrOfLoci);
  stopifnot(length(posOutlierCall) == nbrOfLoci);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Update 'DNAcopy' object
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # (a) segmentation table
  data <- getLocusData(fit);
  data$negOutlierCall <- negOutlierCall;
  data$posOutlierCall <- posOutlierCall;
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



setMethodS3("extractCallsByLocus", "CBS", function(fit, ...) {
  nbrOfLoci <- nbrOfLoci(fit);

  # Extract locus data
  data <- getLocusData(fit);

  # Extract segment data
  segs <- fit$output;

  # Identify segment calls
  callCols <- grep("Call$", colnames(segs));
  nbrOfCalls <- length(callCols);


  chromosome <- data$chromosome;
  x <- data$x;
  y <- data[,3];

  # Allocate locus calls
  naValue <- as.logical(NA);
  callsL <- matrix(naValue, nrow=nbrOfLoci, ncol=nbrOfCalls);
  colnames(callsL) <- colnames(segs)[callCols];
  callsL <- as.data.frame(callsL);

  # For each segment...
  for (ss in seq(length=nrow(segs))) {
    seg <- segs[ss,];
    idxs <- which(chromosome == seg$chromosome & 
                  seg$start <= x & x <= seg$end);
    idxs <- Arguments$getIndices(idxs, max=nbrOfLoci);
    # Sanity check
##    stopifnot(length(idxs) == seg$nbrOfLoci);

    callsSS <- seg[callCols];
    for (cc in seq(length=nbrOfCalls)) {
      callsL[idxs,cc] <- callsSS[,cc];
    }
  } # for (ss ...)

  # The calls for loci that have missing annotations or observations, 
  # should also be missing, i.e. NA.
  nok <- (is.na(chromosome) | is.na(x) | is.na(y));
  callsL[nok,] <- as.logical(NA);

  # Sanity check
  stopifnot(nrow(callsL) == nbrOfLoci);
  stopifnot(ncol(callsL) == nbrOfCalls);

  callsL;
}, private=TRUE) # extractCallsByLocus()



###########################################################################/**
# @RdocMethod getCallStatistics
#
# @title "Calculates various call statistics per chromosome"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{regions}{An optional @data.frame with columns "chromosome", 
#     "start", and "end" specifying the regions of interest to calculate
#     statistics for.  If @NULL, all of the genome is used.}
#  \item{...}{Not used.}
# }
#
# \value{
#  Returns a CxK @data.frame, where C is the number of regions that
#  meet the criteria setup by argument \code{regions} 
#  and (K-4)/2 is the number of call types.
#  The first column is the chromosome index, the second and the third
#  are the first and last position, and the fourth the length
#  (=last-first+1) of the chromosome.
#  The following columns contains call summaries per chromosome.
#  For each chromosome and call type, the total length of such calls
#  on that chromosome is reported together how large of a fraction
#  of the chromosome such calls occupy.
# }
#
# \details{
#   The estimators implemented here are based solely on the 
#   segmentation results, which is very fast.  
#   In the original proposal by Fridlyand et al. [1], the authors
#   estimates the parameters by converting segment-level calls back
#   to locus-level calls and there do the calculations.  
#   The difference between the two approaches should be minor, 
#   particularly for large density arrays.
# }
#
# @author
#
# \references{
#   [1] Fridlyand et al. \emph{Breast tumor copy number aberration 
#       phenotypes and genomic instability}, BMC Cancer, 2006. \cr
# }
#
# \seealso{
#   @seeclass
# }
#
# @keyword internal 
#*/###########################################################################  
setMethodS3("getCallStatistics", "CBS", function(fit, regions=NULL, ...) {
  # To please R CMD check, cf. subset()
  chromosome <- NULL; rm(chromosome);

  # Argument 'regions':
  if (is.null(regions)) {
    # Get chromosome lengths
    regions <- getChromosomeRanges(fit)[,c("chromosome", "start", "end")];
  }
  regions <- as.data.frame(regions);
  stopifnot(all(is.element(c("chromosome", "start", "end"), colnames(regions))));
  stopifnot(!any(duplicated(regions$chromosome)));

  # Calculate lengths
  regions$length <- regions[,"end"] - regions[,"start"] + 1L;

  # Filter out segments within the requested regions
  segsT <- NULL;
  segs <- getSegments(fit, splitter=FALSE);
  for (rr in seq(length=nrow(regions))) {
    regionRR <- regions[rr,];
    chrRR <- regionRR[,"chromosome"];
    startRR <- regionRR[,"start"];
    endRR <- regionRR[,"end"];
    if (is.na(chrRR) || is.na(startRR) || is.na(endRR)) {
      next;
    }

    # Select regions that (at least) overlapping with the region
    segsRR <- subset(segs, chromosome == chrRR & start <= endRR & end >= startRR);

    # Skip?
    if (nrow(segsRR) == 0) {
      next;
    }

    # Adjust ranges
    segsRR$start[segsRR$start < startRR] <- startRR;
    segsRR$end[segsRR$end > endRR] <- endRR;

    segsRR$fullLength <- endRR - startRR + 1L;

    segsT <- rbind(segsT, segsRR);
  } # for (rr ...)

  segs <- segsT;

  # Sum length of calls per type and chromosome
  segs$length <- segs[,"end"] - segs[,"start"] + 1L;
  segs <- segs[order(segs$chromosome),];

  callTypes <- grep("Call$", colnames(segs), value=TRUE);
  res <- lapply(callTypes, FUN=function(type) {
    coeffs <- as.integer(segs[,type]);
    lens <- coeffs * segs$length;
    lens <- by(lens, INDICES=segs$chromosome, FUN=sum, na.rm=TRUE);
    as.vector(lens);
  });
  names(res) <- gsub("Call$", "Length", callTypes);
  res1 <- as.data.frame(res);

  # Extract selected regions
  idxs <- match(unique(segs$chromosome), regions$chromosome);
  regionsT <- regions[idxs,];

  # Sanity check
  stopifnot(nrow(regionsT) == nrow(res1));

  res2 <- res1 / regionsT[,"length"];
  names(res2) <- gsub("Call$", "Fraction", callTypes);

  res3 <- cbind(res1, res2);

  res <- regionsT;
  if (nrow(res3) > 0) {
    res <- cbind(res, res3);
  }
  rownames(res) <- NULL;

  # Sanity checks
  res <- res[,grep("Fraction", colnames(res))];
  for (key in colnames(res)) {
    rho <- res[,key];
    stopifnot(all(rho >= 0, na.rm=TRUE));
    stopifnot(all(rho <= 1, na.rm=TRUE));
  }

  res;
}, protected=TRUE) # getCallStatistics()



###########################################################################/**
# @RdocMethod getFractionOfGenomeLost
# @aliasmethod getFractionOfGenomeGained
# @aliasmethod getFractionOfGenomeAltered
# @aliasmethod getFGL
# @aliasmethod getFGG
# @aliasmethod getFGA
#
# @title "Calculates the fraction of the genome lost, gained, or aberrant either way"
#
# \description{
#  @get "title" (in sense of total copy numbers),
#  using definitions closely related to those presented in [1].
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#  Returns a @double in [0,1].
# }
#
# @author
#
# \references{
#   [1] Fridlyand et al. \emph{Breast tumor copy number aberration 
#       phenotypes and genomic instability}, BMC Cancer, 2006. \cr
# }
#
# \seealso{
#   Internally, @seemethod "getCallStatistics" is used.
#   @seeclass
# }
#
# @keyword internal 
#*/###########################################################################  
setMethodS3("getFractionOfGenomeLost", "CBS", function(fit, ...) {
  stats <- getCallStatistics(fit, ...);
  mean(stats$lossFraction, na.rm=TRUE);
})

setMethodS3("getFractionOfGenomeGained", "CBS", function(fit, ...) {
  stats <- getCallStatistics(fit, ...);
  mean(stats$gainFraction, na.rm=TRUE);
})

setMethodS3("getFractionOfGenomeAltered", "CBS", function(fit, ...) {
  getFractionOfGenomeLost(fit, ...) + getFractionOfGenomeGained(fit, ...);
})

# Shortcuts
setMethodS3("getFGL", "CBS", function(fit, ...) {
  getFractionOfGenomeLost(fit, ...);
}, protected=TRUE)  

setMethodS3("getFGG", "CBS", function(fit, ...) {
  getFractionOfGenomeGained(fit, ...);
}, protected=TRUE)  

setMethodS3("getFGA", "CBS", function(fit, ...) {
  getFractionOfGenomeAltered(fit, ...);
}, protected=TRUE)




setMethodS3("isWholeChromosomeGained", "CBS", function(fit, minFraction=0.99, ...) {
  # Argument 'minFraction':
  minFraction <- Arguments$getDouble(minFraction, range=c(0,1));

  stats <- getCallStatistics(fit, ...);
  calls <- stats$gainFraction;
  if (is.null(calls)) {
    return(rep(NA, times=nbrOfChromosomes(fit)));
  }

  res <- (calls >= minFraction);
  names(res) <- stats$chromosome;
  attr(res, "minFraction") <- minFraction;

  res;
}) # isWholeChromosomeGained()


setMethodS3("isWholeChromosomeLost", "CBS", function(fit, minFraction=0.99, ...) {
  # Argument 'minFraction':
  minFraction <- Arguments$getDouble(minFraction, range=c(0,1));

  stats <- getCallStatistics(fit, ...);
  calls <- stats$lossFraction;
  if (is.null(calls)) {
    return(rep(NA, times=nbrOfChromosomes(fit)));
  }

  res <- (calls >= minFraction);
  names(res) <- stats$chromosome;
  attr(res, "minFraction") <- minFraction;

  res;
}) # isWholeChromosomeLost()


setMethodS3("nbrOfLosses", "CBS", function(fit, ...) {
  stats <- getSegments(fit, ...);
  calls <- stats$lossCall;
  if (is.null(calls)) {
    return(as.integer(NA));
  }
  sum(calls, na.rm=TRUE);
})


setMethodS3("nbrOfGains", "CBS", function(fit, ...) {
  stats <- getSegments(fit, ...);
  calls <- stats$gainCall;
  if (is.null(calls)) {
    return(as.integer(NA));
  }
  sum(calls, na.rm=TRUE);
})


setMethodS3("nbrOfAmplifications", "CBS", function(fit, ...) {
  stats <- getSegments(fit, ...);
  calls <- stats$amplificationCall;
  if (is.null(calls)) {
    return(as.integer(NA));
  }
  sum(calls, na.rm=TRUE);
})



############################################################################
# HISTORY:
# 2011-10-06
# o Added optional argument 'regions' to getCallStatistics() for CBS.
# o Now getCallStatistics() for CBS also returns 'start' and 'end'
#   position of each chromosome.
# 2011-10-03
# o DOCUMENTATION: Added more help pages.
# 2011-10-02
# o DOCUMENTATION: Added an Rdoc help page for getFractionOfGenomeLost(),
#   getFractionOfGenomeGained(), getFractionOfGenomeAltered(), getFGL(), 
#   getFGG() and getFGA().
# 2011-09-05
# o Added getCallStatistics() for CBS.
# 2011-09-04
# o Added extractCallsByLocus() for CBS.
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
