setMethodS3("updateMeansTogether", "CBS", function(fit, idxList, ..., FUN=c("mean", "median"), verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  nbrOfSegments <- nbrOfSegments(fit, splitters=TRUE);

  # Argument 'idxList':
  if (!is.list(idxList)) {
    idxList <- list(idxList);
  }
  idxList <- lapply(idxList, FUN=function(idxs) {
    idxs <- Arguments$getIndices(idxs, max=nbrOfSegments);
    sort(unique(idxs));
  });

  # Argument 'FUN':
  FUN <- match.arg(FUN);
  FUN <- get(FUN, mode="function");

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }
 
  verbose && enter(verbose, "Updating mean level estimates of multiple segments");

  verbose && cat(verbose, "Segments:");
  verbose && str(verbose, idxList);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Extract the data and segmentation results
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  data <- getLocusData(fit);

  segs <- getSegments(fit, splitters=TRUE);

  nbrOfSegments <- nrow(segs);
  verbose && cat(verbose, "Total number of segments: ", nbrOfSegments);

  for (ss in seq(along=idxList)) {
    idxs <- idxList[[ss]];

    fitT <- extractSegments(fit, idxs);
    verbose && cat(verbose, "Number of segments: ", nbrOfSegments(fitT));

    dataT <- getLocusData(fitT);
    segsT <- getSegments(fitT);

    y <- dataT$y;

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # Update the TCN segments
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    verbose && enter(verbose, "Recalculate TCN means");

    # (c) Adjust for missing values
    keep <- which(!is.na(y));
  
    # (d) Update mean
    gamma <- FUN(y[keep]);
 
    # Sanity check
    stopifnot(length(gamma) == 0 || !is.na(gamma));
  
    mus <- c(mean=gamma);
  
    verbose && print(verbose, mus);
    verbose && exit(verbose);

    for (key in names(mus)) {
      segs[idxs,key] <- mus[key];
    }
  } # for (ss ...)

  # Return results
  res <- fit;
  res$output <- segs;

  verbose && exit(verbose);

  res;
}, private=TRUE) # updateMeansTogether()
 


############################################################################
# HISTORY:
# 2011-11-28
# o Added updateMeansTogether() for CBS.
# o Created from PairedPSCBS.updateMeansTogether.R.
############################################################################
