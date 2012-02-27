###########################################################################/**
# @set class=PairedPSCBS
# @RdocMethod callGNL
# @aliasmethod callGainNeutralLoss
#
# @title "Calls segments that are gained, copy neutral, or lost"
#
# \description{
#  @get "title", where copy neutral means having a total copy number
#  that corresponds to the ploidy of the genome.
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
#   Returns a @see "PairedPSCBS" object with added calls.
# }
#
# @author
#
# \seealso{
#   Internally, one of the following methods are used:
#   \code{callGNLByTCNofAB()}.
# }
#
#*/###########################################################################
setMethodS3("callGNL", "PairedPSCBS", function(fit, flavor=c("TCN|AB"), ..., minSize=1, force=FALSE) {
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
    fit <- callGNLByTCNofAB(fit, ..., force=force);
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
}) # callGNL()

setMethodS3("callGainNeutralLoss", "PairedPSCBS", function(...) {
  callGNL(...);
})



setMethodS3("callGNLByTCNofAB", "PairedPSCBS", function(fit, deltaLoss=-0.5, deltaGain=+0.5, alpha=0.05, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'alpha':
  disallow <- c("NA", "NaN", "Inf");
  alpha <- Arguments$getDouble(alpha, range=c(0,0.5), disallow=disallow);

  # Argument 'deltaLoss' & 'deltaGain':
  disallow <- c("NA", "NaN", "Inf");
  deltaLoss <- Arguments$getDouble(deltaLoss, range=c(-Inf,0), disallow=disallow);
  deltaGain <- Arguments$getDouble(deltaGain, range=c(0,+Inf), disallow=disallow);

  # Argument 'force':
  force <- Arguments$getLogical(force);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "callGNLByTCNofAB");
  verbose && cat(verbose, "Alpha: ", alpha);
  verbose && cat(verbose, "Delta loss: ", deltaLoss);
  verbose && cat(verbose, "Delta gain: ", deltaGain);

  segs <- getSegments(fit, splitters=TRUE, simplify=FALSE);

  # Nothing to do?
  if (!force && all(is.element(c("gainCall", "cnCall", "lossCall"), names(segs$cnCall)))) {
    # Segments are already called
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

  verbose && enter(verbose, "Calling segments");


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Confidence interval of copy-neutral AB segments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get call regions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  naValue <- as.double(NA);
  callRegions <- matrix(c(
     Inf,    1,
       1,    1,
       1,  Inf
  ), nrow=3, ncol=2, byrow=TRUE);
  rownames(callRegions) <- c("loss", "cn", "gain");
  colnames(callRegions) <- c("lower", "upper");
  callRegions["loss",] <- ci[1]+callRegions["loss",]*deltaLoss;
  callRegions[  "cn",] <- ci   +callRegions[  "cn",]*c(deltaLoss, deltaGain);
  callRegions["gain",] <- ci[2]+callRegions["gain",]*deltaGain;

  verbose && cat(verbose, "Call (\"acceptance\") regions:");
  verbose && print(verbose, callRegions);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get statistics for all segments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  nbrOfSegs <- nrow(segs);
  verbose && cat(verbose, "Number of segments: ", nbrOfSegs);
  nbrOfABs <- sum(segs$abCall, na.rm=TRUE);
  verbose && cat(verbose, "Number of AB segments: ", nbrOfABs);
  verbose && cat(verbose, "Number of non-AB segments: ", nbrOfSegs-nbrOfABs);

  # Get TCN confidence intervals for all segments
  keys <- sprintf("tcn_%g%%", 100*c(alpha/2, 1-alpha/2));
  ci <- segs[,keys];

  # Call states
  for (rr in seq(length=nrow(callRegions))) {
    state <- rownames(callRegions)[rr];
    verbose && enter(verbose, "Identify all '", state, "' segments");;
    range <- callRegions[rr,];
    verbose && printf(verbose, "Call (\"acceptance\") region: [%g,%g]\n", range[1], range[2]);

    # If a confidence interval is completely within the
    # calling region, call it
    isCalled <- (range[1] <= ci[,1] & ci[,2] < range[2]);

    nbrOfCalled <- sum(isCalled, na.rm=TRUE);
    verbose && cat(verbose, "Number of segments called '", state, "': ", nbrOfCalled);
##    verbose && cat(verbose, "Number of non-AB segments called '", state, "': ", (nbrOfSegs-nbrOfABs)-nbrOfCalled);

    key <- sprintf("%sCall", state);
    segs[[key]] <- isCalled;
    verbose && exit(verbose);
  } # for (rr ...)

  fitC <- fit;
  fitC$output <- segs;

  verbose && exit(verbose);
  
  fitC;
}, protected=TRUE) # callGNLByTCNofAB()



##############################################################################
# HISTORY
# 2012-02-26 [HB]
# o Added internal callGNLByTCNofAB().
# o Added callGNL().
##############################################################################
