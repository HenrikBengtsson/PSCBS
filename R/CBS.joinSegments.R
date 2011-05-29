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
    x <- fit$data$maploc;
    prevSeg <- segs[1L,];
    for (ss in 2:nbrOfSegs) {
      currSeg <- segs[ss,];
      currStart <- currSeg[,"loc.start"];
      prevEnd <- prevSeg[,"loc.end"];
  
      # Center CP
      xMid <- (prevEnd + currStart) / 2;
  
      # Move previous end and current start to this centered CP
      segs[ss,"loc.start"] <- xMid;
      segs[ss-1L,"loc.end"] <- xMid;
  
      prevSeg <- currSeg;
    } # for (ss ...)
    verbose && exit(verbose);
  } # if (nbrOfSegs > 1)


  knownCPs <- range;
  if (!is.null(knownCPs)) {
    if (nbrOfSegs > 0) {
      # Sanity check
      stopifnot(knownCPs[1L] <= segs[1L,"loc.start"]);
      segs[1L,"loc.start"] <- knownCPs[1L];
      # Sanity check
      stopifnot(segs[1L,"loc.end"] <= knownCPs[length(knownCPs)]);
      segs[nbrOfSegs,"loc.end"] <- knownCPs[length(knownCPs)];
    } # if (nbrOfSegs > 0)
  } # if (!is.null(knownCPs))

  fit$output <- segs;

  verbose && print(verbose, head(as.data.frame(fit)));
  verbose && print(verbose, tail(as.data.frame(fit)));
  verbose && exit(verbose);

  fit;
}) # joinSegments()


############################################################################
# HISTORY:
# 2010-11-21
# o Extracted from segmentByPairedPSCBS.R
############################################################################
