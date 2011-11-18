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
  segs <- getSegments(fit);
  verbose && cat(verbose, "Segments:");
  verbose && print(verbose, segs);
  verbose && cat(verbose, "Range:");
  verbose && print(verbose, range);

  nbrOfSegs <- nrow(segs);
  if (nbrOfSegs > 1) {
    verbose && enter(verbose, "Centering change points");
    data <- getLocusData(fit);
    x <- data$x;
    prevSeg <- segs[1L,];
    for (ss in 2:nbrOfSegs) {
      currSeg <- segs[ss,];
      currStart <- currSeg[,"start"];
      prevEnd <- prevSeg[,"end"];

      # Sanity check
      stopifnot(all(currStart >= prevEnd, na.rm=TRUE));
  
      # Center CP
      xMid <- (prevEnd + currStart) / 2;

      # Move previous end and current start to this centered CP
      segs[ss,"start"] <- xMid;
      segs[ss-1L,"end"] <- xMid;

      prevSeg <- currSeg;
    } # for (ss ...)
    verbose && exit(verbose);

    # Sanity checks
    stopifnot(all(segs$start[-1] >= segs$end[-nbrOfSegs], na.rm=TRUE));
    stopifnot(all(diff(segs$start) > 0, na.rm=TRUE));
    stopifnot(all(diff(segs$end) > 0, na.rm=TRUE));

    if (nbrOfSegs > 6) {
      verbose && print(verbose, head(segs));
      verbose && print(verbose, tail(segs));
    } else {
      verbose && print(verbose, segs);
    }
  } # if (nbrOfSegs > 1)


  if (!is.null(range)) {
    verbose && enter(verbose, "Adjust for 'range'");
    verbose && cat(verbose, "Range:");
    verbose && print(verbose, range);
    xMin <- min(range, na.rm=TRUE);
    xMax <- max(range, na.rm=TRUE);
    if (nbrOfSegs > 0) {
      # Sanity check
      stopifnot(xMin <= segs[1L,"start"]);
      segs[1L,"start"] <- xMin;
      # Sanity check
      stopifnot(segs[1L,"end"] <= xMax);
      segs[nbrOfSegs,"end"] <- xMax;

      # Sanity checks
      stopifnot(all(segs$start[-1] >= segs$end[-nbrOfSegs], na.rm=TRUE));
      stopifnot(all(diff(segs$start) > 0, na.rm=TRUE));
      stopifnot(all(diff(segs$end) > 0, na.rm=TRUE));


      if (nbrOfSegs > 6) {
        verbose && print(verbose, head(segs));
        verbose && print(verbose, tail(segs));
      } else {
        verbose && print(verbose, segs);
      }
    } # if (nbrOfSegs > 0)
    verbose && exit(verbose);
  } # if (!is.null(range))

  fit$output <- segs;

  segs <- as.data.frame(fit);
  if (nbrOfSegs > 6) {
    verbose && print(verbose, head(segs));
    verbose && print(verbose, tail(segs));
  } else {
    verbose && print(verbose, segs);
  }
  verbose && exit(verbose);

  fit;
}, private=TRUE) # joinSegments()


############################################################################
# HISTORY:
# 2011-11-17
# o Added more sanity checks to joinSegments().
# 2011-09-04
# o Updated joinSegments() to be aware of new column names in CBS.
# 2011-06-14
# o Updated code to recognize new column names.
# 2010-11-21
# o Extracted from segmentByPairedPSCBS.R
############################################################################
