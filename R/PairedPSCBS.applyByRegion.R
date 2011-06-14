setMethodS3("applyByRegion", "PairedPSCBS", function(fit, FUN, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'FUN':
  stopifnot(is.function(FUN));

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  } 


  verbose && enter(verbose, "Apply function region by region");


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract data and estimates
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  data <- fit$data;
  tcnSegRows <- fit$tcnSegRows;
  dhSegRows <- fit$dhSegRows;
  segs <- fit$output;
  params <- fit$params;

  # Sanity checks
  if (!params$joinSegments) {
    throw("Cannot applyByRegion() unless PSCNs are segmented using joinSegments=TRUE.");
  }
  dataRows <- tcnSegRows;

  # Sanity checks
  stopifnot(all(!is.na(data$chromosome) & !is.na(data$x)));
  stopifnot(length(tcnSegRows) == length(dhSegRows));


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # For each segment...
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  nbrOfSegments <- nrow(segs);
  verbose && cat(verbose, "Number of segments: ", nbrOfSegments);

  # Allocate result objects
  dataN <- outputN <- dataRowsN <- NULL;

  for (rr in seq(length=nbrOfSegments)) {
    verbose && enter(verbose, sprintf("Segment #%d of %d", rr, nbrOfSegments));

    # Extract segment
    segRR <- segs[rr,,drop=FALSE];

    # Nothing todo?
    if (is.na(segRR[["tcnId"]]) && is.na(segRR[["dhId"]])) {
      verbose && cat(verbose, "A divider. Nothing to do.");
      outputN <- rbind(outputN, NA);
      dataRowsN <- rbind(dataRowsN, NA);
      verbose && exit(verbose);
      next;
    }

    verbose && str(verbose, segRR, level=-20);

    # Extract data
    dataRowsRR <- dataRows[rr,,drop=FALSE];
    from <- dataRowsRR[[1]];
    to <- dataRowsRR[[2]];
    ok <- (!is.na(from) & !is.na(to));
    from <- from[ok];
    to <- to[ok];
    keep <- logical(nrow(data));
    for (kk in seq(along=from)) {
      keep[from[kk]:to[kk]] <- TRUE;
    }
    dataRowsRR <- which(keep);
    verbose && printf(verbose, "Identified %d (%.2f%%) of %d data rows:\n", length(dataRowsRR), 100*length(dataRowsRR)/nrow(data), nrow(data));
    verbose && str(verbose, dataRowsRR);
    dataRR <- data[dataRowsRR,,drop=FALSE];
    verbose && str(verbose, dataRR, level=-20);

    verbose && enter(verbose, "Applying function 'FUN' to segment");
    resRR <- FUN(kk, segRR, dataRR, ...);
    verbose && cat(verbose, "Returned result:");
    verbose && str(verbose, resRR, level=-20);
    verbose && exit(verbose);

    # Nothing to update/store?
    if (is.null(resRR)) {
      verbose && cat(verbose, "Segment and its data was dropped: ", rr);
      verbose && exit(verbose);
      next;
    }

    # Update locus-level data?
    dataRRN <- resRR$data;
    dataRowsRRN <- c(1L, nrow(dataRRN));
    if (!is.null(dataN)) {
      dataRowsRRN <- dataRowsRRN + nrow(dataN);
    }
    dataN <- rbind(dataN, dataRRN);

    # Update segment table?
    segRRN <- resRR$output;
    outputN <- rbind(outputN, segRRN);
    dataRowsN <- rbind(dataRowsN, dataRowsRRN);

    # Sanity checks
    stopifnot(nrow(dataN) == max(dataRowsN, na.rm=TRUE));
    stopifnot(nrow(outputN) == nrow(dataRowsN));

    verbose && exit(verbose);
  } # for (rr ...)

  rownames(dataRowsN) <- NULL;
  colnames(dataRowsN) <- colnames(dataRows);
  dataRowsN <- as.data.frame(dataRowsN);

  # Sanity checks
  stopifnot(!is.null(dataN));
  stopifnot(!is.null(outputN));
  stopifnot(!is.null(dataRowsN));

  # Create result
  res <- fit; # "clone"
  res$data <- dataN;
  res$output <- outputN;
  res$tcnSegRows <- dataRowsN;
  res$dhSegRows <- dataRowsN;  # Is this really the case? /HB 2011-01-17

  verbose && exit(verbose);

  res;
}, private=TRUE)


.addC1C2WithStatitics <- function(rr, output, data, robust=TRUE, ...) {
  # Calculate locus-level (C1,C2)
  C <- data$CT;
  rho <- data$rho;
  C1 <- 1/2 * (1 - rho) * C;
  C2 <- C - C1;
  CC <- data.frame(C1=C1, C2=C2);

  if (robust) {
    meanFcn <- function(x, ...) median(x, na.rm=TRUE);
    sdFcn <- function(x, ...) mad(x, na.rm=TRUE);
  } else {
    meanFcn <- function(x, ...) mean(x, na.rm=TRUE);
    sdFcn <- function(x, ...) sd(x, na.rm=TRUE);
  }

  # Calculate region-level (C1,C2) means and std devs.
  muCC <- apply(CC, MARGIN=2, FUN=meanFcn);
  sigmaCC <- apply(CC, MARGIN=2, FUN=sdFcn);
  rhoCC <- cor(CC[,1], CC[,2], use="pairwise.complete.obs");

  names(muCC) <- c("c1Avg", "c2Avg");
  names(sigmaCC) <- c("c1Sd", "c2Sd");

  # Update data
  data <- cbind(data, CC);


  # Update segment table
  outputT <- c(muCC, sigmaCC, c1c2.cor=rhoCC);
  outputT <- as.list(outputT);
  outputT <- as.data.frame(outputT);

  output <- cbind(output, outputT);

  list(data=data, output=output);
} # .addC1C2WithStatitics()


.addCACBWithStatitics <- function(rr, output, data, beta=c("betaTN", "betaT"), stratifyBy=c("all", "hets", "homs"), robust=TRUE, ...) {
  # Argument 'beta':
  beta <- match.arg(beta);

  # Argument 'stratifyBy':
  stratifyBy <- match.arg(stratifyBy);

  # Calculate locus-level (CA,CB)
  C <- data$CT;
  beta <- data[[beta]];
  CB <- beta * C;
  CA <- C - CB;
  CC <- data.frame(CA=CA, CB=CB);

  # Update data
  data <- cbind(data, CC);


  if (robust) {
    meanFcn <- function(x, ...) median(x, na.rm=TRUE);
    sdFcn <- function(x, ...) mad(x, na.rm=TRUE);
  } else {
    meanFcn <- function(x, ...) mean(x, na.rm=TRUE);
    sdFcn <- function(x, ...) sd(x, na.rm=TRUE);
  }

  if (stratifyBy == "hets") {
    muN <- data$muN;
    keep <- (muN == 1/2);
    CC <- CC[keep,,drop=FALSE];
  } else if (stratifyBy == "homs") {
    muN <- data$muN;
    keep <- (muN == 0 | muN == 1);
    CC <- CC[keep,,drop=FALSE];
  }

  # Calculate region-level (CA,CB) means and std devs.
  muCC <- apply(CC, MARGIN=2, FUN=meanFcn);
  sigmaCC <- apply(CC, MARGIN=2, FUN=sdFcn);
  if (nrow(CC) < 3) {
    rhoCC <- as.double(NA);
  } else {
    rhoCC <- cor(CC[,1], CC[,2], use="pairwise.complete.obs");
  }
  names(muCC) <- c("caAvg", "cbAvg");
  names(sigmaCC) <- c("caSd", "cbSd");

  # Update segment table
  outputT <- c(muCC, sigmaCC, cacbCor=rhoCC);
  outputT <- as.list(outputT);
  outputT <- as.data.frame(outputT);

  output <- cbind(output, outputT);

  list(data=data, output=output);
} # .addCACBWithStatitics()



#############################################################################
# HISTORY:
# 2011-06-14
# o Updated code to recognize new column names.
# 2011-01-27
# o Added .addCACBWithStatitics().
# o Added .addC1C2WithStatitics().
# o Added applyByRegion().
# o Created.
#############################################################################
