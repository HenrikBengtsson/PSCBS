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

  segs <- getSegments(fit, splitters=TRUE);

  # Nothing to do?
  if (!force && !is.null(segs$cnCall)) {
    # Copy neutral segments are already called
    return(fit);
  }

  # Getting AB calls
  isAB <- segs$abCall;
  if (is.null(isAB)) {
    throw("Cannot call copy-neutral states, because allelic-balance calls have not been made yet.");
  }

  C <- segs[,"tcnMean", drop=TRUE];
  isAB <- segs[,"abCall", drop=TRUE];
  n <- segs[,"tcnNbrOfSNPs", drop=TRUE]; # "tcnNbrOfLoci"? /HB 2010-09-09

  # Give more weight to longer regions
  weights <- n;

  isNeutral <- findNeutralCopyNumberState(C=C, isAI=!isAB, weights=weights,
                                                       ..., verbose=verbose);

  segs$cnCall <- isNeutral;

  fitC <- fit;
  fitC$output <- segs;
  
  fitC;
}, protected=TRUE) # callCopyNeutralByTCNofAB()




##############################################################################
# HISTORY
# 2012-02-24 [HB]
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
