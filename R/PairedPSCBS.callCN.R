###########################################################################/**
# @set class=PairedPSCBS
# @RdocMethod callCN
# @aliasmethod callCopyNeutral
#
# @title "Calls segments that are copy neutral"
#
# \description{
#  @get "title", i.e. that have a total copy number that corresponds
#  to the ploidy of the genome.
# }
#
# @synopsis
#
# \arguments{
#   \item{flavor}{A @character string specifying which type of
#    call to use.}
#   \item{...}{Additional arguments passed to the caller.}
#   \item{minSize}{An optional @integer specifying the minimum number
#    of data points in order to call a segments.  If fewer data points,
#    then the call is set to @NA regardless.}
#   \item{force}{If @FALSE, and copy-neutral calls already exits,
#    then nothing is done, otherwise the calls are done.}
# }
#
# \value{
#   Returns a @see "PairedPSCBS" object with copy-neutral calls.
# }
#
# @author
#
# \seealso{
#   Internally, one of the following methods are used:
#   @seemethod "callCopyNeutralByTCNofAB".
# }
#
#*/###########################################################################
setMethodS3("callCN", "PairedPSCBS", function(fit, flavor=c("TCN|AB"), ..., minSize=1, force=FALSE) {
  # Argument 'flavor':
  flavor <- match.arg(flavor);

  # Argument 'minSize':
  minSize <- Arguments$getDouble(minSize, range=c(1,Inf));


  # Already done?
  segs <- as.data.frame(fit);
  calls <- segs$cnCall;
  if (!force && !is.null(calls)) {
    return(invisible(fit));
  }
  
  if (flavor == "TCN|AB") {
    fit <- callCopyNeutralByTCNofAB(fit, ...);
  } else {
    throw("Cannot call allelic balance. Unsupported flavor: ", flavor);
  }

  # Don't call segments with too few data points?
  if (minSize > 1) {
    segs <- as.data.frame(fit);
    ns <- segs$dhNbrOfLoci;
    calls <- segs$cnCall;
    calls[ns < minSize] <- NA;
    segs$cnCall <- calls;
    fit$output <- segs;
    rm(segs, ns, calls); # Not needed anymore
  }

  return(invisible(fit));
})

setMethodS3("callCopyNeutral", "PairedPSCBS", function(...) {
  callCN(...);
})



###########################################################################/**
# @set class=PairedPSCBS
# @RdocMethod callCopyNeutralByTCNofAB
#
# @title "Calls regions that are copy neutral"
#
# \description{
#  @get "title" from the total copy numbers (TCNs) of segments
#  in allelic balance (AB).
# }
#
# @synopsis
#
# \arguments{
#   \item{fit}{A PairedPSCBS fit object as returned by 
#     @see "PSCBS::segmentByPairedPSCBS".}
#   \item{...}{Additional arguments passed to 
#     @see "findNeutralCopyNumberState".}
#   \item{force}{If @TRUE, an already called object is skipped, otherwise not.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   Returns a @see "PairedPSCBS" fit object where a column 
#   with the copy-neutral call.
# }
#
# %% examples "../incl/callCopyNeutralByTCNofAB.PairedPSCBS.Rex"
#
# @author
#
# @keyword internal
#*/########################################################################### 
setMethodS3("callCopyNeutralByTCNofAB", "PairedPSCBS", function(fit, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "callCopyNeutralByTCNofAB");

  segs <- getSegments(fit, splitters=TRUE);

  # Nothing to do?
  if (!force && !is.null(segs$cnCall)) {
    # Copy neutral segments are already called
    return(fit);
  }


  verbose && enter(verbose, "Identifying copy neutral AB segments");

  # Getting AB calls
  isAB <- segs$abCall;
  if (is.null(isAB)) {
    throw("Cannot call copy-neutral states, because allelic-balance calls have not been made yet.");
  }

  nABs <- sum(isAB, na.rm=TRUE);
  verbose && cat(verbose, "Number of AB segments: ", nABs);
  if (nABs == 0) {
    throw("Cannot call copy-neutral states, because none of the segments are in allelic balance.");
  }

  C <- segs[,"tcnMean", drop=TRUE];
  isAB <- segs[,"abCall", drop=TRUE];
  n <- segs[,"tcnNbrOfSNPs", drop=TRUE]; # "tcnNbrOfLoci"? /HB 2010-09-09

  # Give more weight to longer regions
  weights <- n;

  # Identify copy neutral AB segments
  isNeutralAB <- findNeutralCopyNumberState(C=C, isAI=!isAB, weights=weights,
                                                       ..., verbose=verbose);
  n <- sum(isNeutralAB, na.rm=TRUE);
  verbose && cat(verbose, "Number of copy-neutral AB segments: ", n);
  if (n == 0) {
    throw("Cannot call copy-neutral states, because none of the segments in allelic-balance are copy neutral.");
  }

  verbose && enter(verbose, "Extracting all copy neutral AB segments across all chromosomes into one big segments");

  # Estimate common quantiles for those segments
  fitCN <- extractSegments(fit, isNeutralAB);
  # Fake one chromosome
  fitCN$data$chromosome <- 0L;
  fitCN$output$chromosome <- 0L;
  verbose && print(verbose, fitCN);
  verbose && exit(verbose);

  # Turn into one big segment by dropping all change points
  nCPs <- nbrOfChangePoints(fitCN, splitters=TRUE);
  if (nCPs > 0) {
    for (kk in nCPs:1) {
      fitCN <- dropChangePoint(fitCN, kk, update=FALSE);
    } # for (kk ...)
    fitCN <- updateMeans(fitCN);
  }
  # Sanity check
  stopifnot(nbrOfSegments(fitCN) == 1);
  verbose && exit(verbose);

  verbose && enter(verbose, "Estimating bootstrap TCN quantiles for the copy-neutral state");
  fitCN <- bootstrapTCNandDHByRegion(fitCN, force=TRUE, ..., 
                                     verbose=less(verbose, 10));
  segsCN <- getSegments(fitCN, simplified=FALSE);
  keep <- grep("^tcn_.*%$", colnames(segsCN));
  segsCN <- segsCN[,keep,drop=FALSE];
  verbose && print(verbose, segsCN);
  # Sanity check
  stopifnot(ncol(segsCN) > 0);
  verbose && exit(verbose);

  # Call all segments
  alpha <- 0.05/2;
  keys <- sprintf("tcn_%g%%", 100*c(alpha, 1-alpha));
  segsCN <- segsCN[,keys];
  # Sanity check
  stopifnot(ncol(segsCN) > 0);
  range <- unlist(segsCN[1,]);
  verbose && print(verbose, range);

  tcnMean <- segs$tcnMean;
  isNeutral <- (range[1] <= tcnMean & tcnMean <= range[2]);

  n <- sum(isNeutral, na.rm=TRUE);
  verbose && cat(verbose, "Number of segments called copy-neutral: ", n);
  
  # Sanity check
  # All previously called AB regions should remain called here as well
  stopifnot(all(isNeutral[isNeutralAB], na.rm=TRUE));

  segs$cnCall <- isNeutral;

  fitC <- fit;
  fitC$output <- segs;

  verbose && exit(verbose);
  
  fitC;
}, protected=TRUE) # callCopyNeutralByTCNofAB()




##############################################################################
# HISTORY
# 2012-02-24 [HB]
# o Now callCopyNeutralByTCNofAB() calls all segements, not just those in AB.
# o Now the copy-neutral calls are named 'cnCall' (not 'neutralCall').
# o Added callCN()/callCopyNeutral().
# o Added callCopyNeutralByTCNofAB() for PairedPSCBS.  The method was 
#   adopted from callCopyNeutralRegions() in aroma.cn, whose history has 
#   been incorporated below.
# o Created.
# 2010-09-15* [HB]
# o Added Rdocs for callCopyNeutralRegions(). 
# 2010-09-09* [HB]
# o Added callCopyNeutralRegions() for PairedPSCBS.
##############################################################################
