##########################################################################/**
# @set class=PairedPSCBS
# @RdocMethod callROH
#
# @title "Calls segments that are in ROH"
#
# \description{
#  @get "title", i.e. that have no (true) heterozygous genotypes.
#  Run of homozygosity (ROH) is a property of the normal (germline) sample.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Additional arguments passed to @see "testROH".}
#   \item{force}{If @FALSE, and ROH calls already exits,
#    then nothing is done, otherwise the calls are done.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   Returns a @see "PairedPSCBS" object with ROH calls.
# }
#
# @author
#
# \seealso{
#   Internally, @see "testROH" is used.
#   To call allelic balance (AB) see @seemethod "callAB".
#   To call loss of heterozygosity (LOH) see @seemethod "callLOH".
# }
#*/########################################################################### 
setMethodS3("callROH", "PairedPSCBS", function(fit, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Calling ROH");

  # Already done?
  segs <- getSegments(fit);
  calls <- segs$rohCall;
  if (!force && !is.null(calls)) {
    return(invisible(fit));
  }

  nbrOfSegments <- nrow(segs);
  calls <- rep(NA, times=nbrOfSegments);
  if (is.null(calls)) {
    segs <- cbind(segs, rohCall=calls);
  }

  # For each segment...
  for (ss in seq(length=nbrOfSegments)) {
    verbose && enter(verbose, sprintf("Segment #%d of %d", ss, nbrOfSegments));

    fitT <- extractSegment(fit, ss);

    # Call only "non-splitter" segments
    if (nbrOfSegments(fitT) > 0) {
      calls[ss] <- callROHOneSegment(fitT, ..., verbose=less(verbose, 1));
    }

    verbose && exit(verbose);
  } # for (ss ...)

  verbose && cat(verbose, "ROH calls:");
  verbose && str(verbose, calls);
  verbose && print(verbose, summary(calls));

  segs$rohCall <- calls;

  fit$output <- segs;

  # Append parameters
  params <- fit$params;
  # To add...
  fit$params <- params;
 
  verbose && exit(verbose);

  invisible(fit);
})




# This method calls ROH for a single-segment PairedPSCBS object
setMethodS3("callROHOneSegment", "PairedPSCBS", function(fit, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Calling ROH for a single segment");

  # Make sure that there is only a single segment in this object
  stopifnot(nbrOfSegments(fit, splitters=TRUE) == 1);


  # Extract the locus-level data for the segment tested
  data <- getLocusData(fit);

  # Keep only SNPs:
  # SNPs are identifies as those loci that have non-missing
  # 'betaTN' & 'muN', cf. segmentByPairedPSCBS().
  isSnp <- (!is.na(data$betaTN) & !is.na(data$muN));
  nbrOfSnps <- sum(isSnp);
  verbose && cat(verbose, "Number of SNPs: ", nbrOfSnps);
  data <- data[isSnp,];

  # Extract that SNP signals used for calling ROH
  betaN <- data$betaN;
  muN <- data$muN;
  csN <- data$csN;  # Genotyping confidence scores, if available

  # Test for ROH
  fit <- testROH(betaN=betaN, muN=muN, csN=csN, ..., verbose=less(verbose, 10));

  # Get the ROH call (TRUE, FALSE, or NA)
  call <- fit;

  verbose && exit(verbose);

  call;
}, private=TRUE)


##############################################################################
# HISTORY
# 2011-11-12 [HB]
# o BUG FIX: ROH calls should be stored in column 'rohCall' (not 'rohCalls').
# 2011-11-04 [HB]
# o Added callROH() for PairedPSCBS.
# o Created.
##############################################################################
