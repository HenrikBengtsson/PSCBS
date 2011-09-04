setMethodS3("joinSegments", "CBS", function(fit, range=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'range':
  if (!is.null(range)) {
    range <- Arguments$getDoubles(range, length=c(2,2));
    stopifnot(range[2] >= range[1]);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Joining segments");  
  segs <- fit$output;

  nbrOfSegs <- nrow(segs);
  if (nbrOfSegs > 1) {
    verbose && enter(verbose, "Centering change points");
    x <- fit$data$x;
    prevSeg <- segs[1L,];
    for (ss in 2:nbrOfSegs) {
      currSeg <- segs[ss,];
      currStart <- currSeg[,"start"];
      prevEnd <- prevSeg[,"end"];
  
      # Center CP
      xMid <- (prevEnd + currStart) / 2;
  
      # Move previous end and current start to this centered CP
      segs[ss,"start"] <- xMid;
      segs[ss-1L,"end"] <- xMid;
  
      prevSeg <- currSeg;
    } # for (ss ...)
    verbose && exit(verbose);
  } # if (nbrOfSegs > 1)


  knownCPs <- range;
  if (!is.null(knownCPs)) {
    if (nbrOfSegs > 0) {
      # Sanity check
      stopifnot(knownCPs[1L] <= segs[1L,"start"]);
      segs[1L,"start"] <- knownCPs[1L];
      # Sanity check
      stopifnot(segs[1L,"end"] <= knownCPs[length(knownCPs)]);
      segs[nbrOfSegs,"end"] <- knownCPs[length(knownCPs)];
    } # if (nbrOfSegs > 0)
  } # if (!is.null(knownCPs))

  fit$output <- segs;

  verbose && print(verbose, head(as.data.frame(fit)));
  verbose && print(verbose, tail(as.data.frame(fit)));
  verbose && exit(verbose);

  fit;
}, private=TRUE) # joinSegments()


############################################################################
# HISTORY:
# 2011-09-04
# o Updated joinSegments() to be aware of new column names in CBS.
# 2011-06-14
# o Updated code to recognize new column names.
# 2010-11-21
# o Extracted from segmentByPairedPSCBS.R
############################################################################
