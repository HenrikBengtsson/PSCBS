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
    fit <- callCopyNeutralByTCNofAB(fit, ..., force=force);
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




setMethodS3("calcStatsForCopyNeutralABs", "PairedPSCBS", function(fit, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'force':
  force <- Arguments$getLogical(force);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "calcStatsForCopyNeutralABs");

  segsCN <- fit$params$copyNeutralStats;
  if (!force && !is.null(segsCN)) {
    verbose && exit(verbose);
    return(fit);
  }

  verbose && enter(verbose, "Identifying copy neutral AB segments");

  # Getting AB calls
  segs <- getSegments(fit, splitters=TRUE);
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
  nAB <- sum(isNeutralAB, na.rm=TRUE);
  verbose && cat(verbose, "Number of copy-neutral AB segments: ", nAB);
  if (nAB == 0) {
    throw("Cannot call copy-neutral states, because none of the segments in allelic-balance are copy neutral.");
  }

  verbose && enter(verbose, "Extracting all copy neutral AB segments across all chromosomes into one big segments");

  # (a) Extract those
  fitCN <- extractSegments(fit, isNeutralAB);
  verbose && print(verbose, fitCN);
  verbose && exit(verbose);

  # (b) Turn into a single-chromosome data set
  fitCN <- extractSegments(fitCN, !isSegmentSplitter(fitCN));

  fitCN$output$chromosome <- 0L;
  fitCN$data$chromosome <- 0L;

  # (c) Turn into one big segment by dropping all change points
  nCPs <- nbrOfChangePoints(fitCN);
  if (nCPs > 1) {
    verbose && enter(verbose, "Dropping all change points");
    fitCN <- dropChangePoints(fitCN, idxs=nCPs:1, update=TRUE, 
                              verbose=less(verbose, 5));
    verbose && exit(verbose);
  }
  # Sanity check
  stopifnot(nbrOfSegments(fitCN) == 1);
  verbose && exit(verbose);

  verbose && enter(verbose, "Bootstrap the identified copy-neutral states");
  fitCN <- bootstrapTCNandDHByRegion(fitCN, force=TRUE, ..., 
                                     verbose=less(verbose, 10));
  segsCN <- getSegments(fitCN, simplified=FALSE);
  names <- colnames(segsCN);
  excl <- grep("(^chromosome|Id|Start|End|Call)$", names);
  segsCN <- segsCN[,-excl,drop=FALSE];
  # Sanity check
  stopifnot(ncol(segsCN) > 0);
  verbose && exit(verbose);

  verbose && print(verbose, segsCN);
  verbose && exit(verbose);

  fit$params$copyNeutralStats <- segsCN;

  invisible(fit);
}, protected=TRUE) # calcStatsForCopyNeutralABs()



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
#   \item{delta}{A non-negative @double specifying the width of the 
#     "acceptance" region.}
#   \item{alpha}{A @double in [0,0.5] specifying the significance level
#     of the confidence intervals used.}
#   \item{...}{Additional arguments passed to 
#              @seemethod "calcStatsForCopyNeutralABs".}
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
setMethodS3("callCopyNeutralByTCNofAB", "PairedPSCBS", function(fit, delta=0.5, alpha=0.05, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'alpha':
  disallow <- c("NA", "NaN", "Inf");
  alpha <- Arguments$getDouble(alpha, range=c(0,0.5), disallow=disallow);

  # Argument 'delta':
  disallow <- c("NA", "NaN", "Inf");
  delta <- Arguments$getDouble(delta, range=c(0,Inf), disallow=disallow);

  # Argument 'force':
  force <- Arguments$getLogical(force);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "callCopyNeutralByTCNofAB");
  verbose && cat(verbose, "Alpha: ", alpha);
  verbose && cat(verbose, "Delta TCN: ", delta);

  segs <- getSegments(fit, splitters=TRUE, simplify=FALSE);

  # Nothing to do?
  if (!force && !is.null(segs$cnCall)) {
    # Copy neutral segments are already called
    verbose && cat(verbose, "Already called. Skipping.");
    verbose && exit(verbose);
    return(fit);
  }

  # Check that bootstrap estimates exists
  keys <- sprintf("tcn_%g%%", 100*c(alpha/2, 1-alpha/2));
  missing <- keys[!is.element(keys, colnames(segs))];
  if (length(missing) > 0) {
    throw("No such statistics: ", hpaste(missing));
  }

  verbose && enter(verbose, "Calling copy-neutral segments");

  verbose && enter(verbose, "Estimating TCN confidence interval of copy-neutral AB segments");

  fit <- calcStatsForCopyNeutralABs(fit, ..., verbose=less(verbose, 5));
  stats <- fit$params$copyNeutralStats;
  verbose && cat(verbose, "Bootstrap statistics for copy-neutral AB segments:");
  verbose && print(verbose, stats);

  # Extract TCN stats
  cols <- grep("^(tcn_|tcnMean)", colnames(stats));
  tcnStats <- stats[,cols,drop=FALSE];
  tcnStats <- unlist(tcnStats, use.names=TRUE);
  verbose && print(verbose, "TCN statistics:");
  verbose && print(verbose, tcnStats);

  # Extract confidence interval of interest
  keys <- sprintf("tcn_%g%%", 100*c(alpha/2, 1-alpha/2));
  missing <- keys[!is.element(keys, names(tcnStats))];
  if (length(missing) > 0) {
    throw("No such statistics: ", hpaste(missing));
  }
  mean <- tcnStats["tcnMean"];
  ci <- tcnStats[keys];
  verbose && printf(verbose, "%g%%-confidence interval of TCN mean for the copy-neutral state: [%g,%g] (mean=%g)\n", 100*(1-alpha), ci[1], ci[2], mean);

  verbose && exit(verbose);

  verbose && enter(verbose, "Identify all copy-neutral segments");;
  range <- ci + delta*c(-1,+1);
  verbose && printf(verbose, "Call (\"acceptance\") region: [%g,%g]\n", range[1], range[2]);

  # Get TCN confidence intervals for all segments
  keys <- sprintf("tcn_%g%%", 100*c(alpha/2, 1-alpha/2));
  ci <- segs[,keys];
  ci <- as.matrix(ci);

  # If a confidence interval is completely within the
  # calling region, call it
  isCN <- (range[1] <= ci[,1] & ci[,2] <= range[2]);

  nbrOfSegs <- nrow(segs);
  nbrOfABs <- sum(segs$abCall, na.rm=TRUE);
  nbrOfCNs <- sum(isCN, na.rm=TRUE);
  verbose && cat(verbose, "Total number of segments: ", nbrOfSegs);
  verbose && cat(verbose, "Number of segments called copy neutral: ", nbrOfCNs);
  verbose && cat(verbose, "Number of non-AB segments called copy neutral: ", (nbrOfSegs-nbrOfABs)-nbrOfCNs);

  verbose && exit(verbose);
  
  # Sanity check
#  # All previously called AB regions should remain called here as well
#  stopifnot(all(isCN[isNeutralAB], na.rm=TRUE));

  segs$cnCall <- isCN;

  fitC <- fit;
  fitC$output <- segs;

  verbose && exit(verbose);
  
  fitC;
}, protected=TRUE) # callCopyNeutralByTCNofAB()



##############################################################################
# HISTORY
# 2012-02-25 [HB]
# o Added internal calcStatsForCopyNeutralABs() for PairedPSCBS.
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
