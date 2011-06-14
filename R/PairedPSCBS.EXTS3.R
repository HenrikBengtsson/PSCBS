setMethodS3("mergeTwoSegments", "PairedPSCBS", function(this, left, ...) {
  nbrOfSegments <- nbrOfSegments(this);
  # Argument 'left':
  left <- Arguments$getIndex(left, max=nbrOfSegments-1L);

  segs <- this$output;
  tcnSegRows <- this$tcnSegRows;
  dhSegRows <- this$dhSegRows;

  rows <- c(left,left+1);
  segsT <- segs[rows,,drop=FALSE];

  # Sanity check
  chrs <- segsT[["chromosome"]];
  stopifnot(chrs[1] == chrs[2]);

  # Merge segments
  segT <- segsT[1,];
  fields <- colnames(segsT);

  # (chromosome, tcnId, dhId)
  idxsUsed <- 1:3;

  # Starts
  idxs <- grep("Start$", fields);
  segT[,idxs] <- apply(segsT[,idxs], MARGIN=2, FUN=min, na.rm=TRUE);
  idxsUsed <- c(idxsUsed, idxs);

  # Ends
  idxs <- grep("End$", fields);
  segT[,idxs] <- apply(segsT[,idxs], MARGIN=2, FUN=max, na.rm=TRUE);
  idxsUsed <- c(idxsUsed, idxs);

  # Counts
  idxs <- grep("NbrOf", fields);
  segT[,idxs] <- apply(segsT[,idxs], MARGIN=2, FUN=sum);
  idxsUsed <- c(idxsUsed, idxs);

  # "Invalidate" remaining entries
  idxsTodo <- setdiff(seq(along=fields), idxsUsed);
  segT[,idxsTodo] <- NA;

  # Update segment table
  segs[rows[1],] <- segT;
  segs <- segs[-rows[2],];

  # Update 'segRows' tables
  segRows <- tcnSegRows;
  segRows[rows[1],2] <- segRows[rows[2],2];
  segRows <- segRows[-rows[2],];
  tcnSegRows <- segRows;

  segRows <- dhSegRows;
  segRows[rows[1],2] <- segRows[rows[2],2];
  segRows <- segRows[-rows[2],];
  dhSegRows <- segRows;

  # Create results object
  res <- this;
  res$output <- segs;
  res$tcnSegRows <- tcnSegRows;
  res$dhSegRows <- dhSegRows;

  # Update the mean estimates.
  res <- updateMeans(res);

  res;
}, private=TRUE)


setMethodS3("dropByRegions", "PairedPSCBS", function(this, regions, H=1, ...) {
  nbrOfSegments <- nbrOfSegments(this);
  # Argument 'regions':
  regions <- Arguments$getIndices(regions, max=nbrOfSegments);

  # Argument 'H':
  H <- Arguments$getInteger(H, range=c(0,Inf));

  # Nothing to do?
  if (H == 0) {
    return(this);
  } else if (H > 1) {
    # Identify regions to drop
    Hs <- seq(length=H);
    regions <- regions - 1L;
    regions <- as.list(regions);
    regions <- lapply(regions, FUN=function(region) region + Hs);
    regions <- unlist(regions, use.names=FALSE);
    regions <- unique(regions);  
  }

  # Identify regions to keep
  allRegions <- seq(length=nbrOfSegments);
  keepRegions <- setdiff(allRegions, regions);

  res <- extractByRegions(this, regions=keepRegions, ...);
  res$dropped <- extractByRegions(this, regions=regions, ...);

  res;
}, private=TRUE)




setMethodS3("extractLocusLevelC1C2", "PairedPSCBS", function(fit, ...) {
  data <- fit$data;
  C <- data$CT;
  rho <- data$rho;
  C1 <- 1/2*(1-rho)*C;
  C2 <- C-C1;
  data.frame(C1=C1, C2=C2);
}, private=TRUE) # extractLocusLevelC1C2()


setMethodS3("extractLocusLevelTCN", "PairedPSCBS", function(fit, ...) {
  data <- fit$data;
  C <- data$CT;
}, private=TRUE) # extractLocusLevelTCN()


setMethodS3("updateMeans", "PairedPSCBS", function(fit, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }
 
  verbose && enter(verbose, "Updating mean level estimates");


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Extract the data and segmentation results
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  data <- fit$data;

  segs <- fit$output;
  keep <- is.finite(segs$chromosome);
  segs <- segs[keep,,drop=FALSE];

  nbrOfSegments <- nrow(segs);
  verbose && cat(verbose, "Number of segments: ", nbrOfSegments);

  chromosome <- data$chromosome;
  x <- data$x;
  CT <- data$CT;
  rho <- data$rho;
  muN <- data$muN;
  isSnp <- !is.na(muN);
  isHet <- isSnp & (muN == 1/2);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Update the TCN segments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  for (ss in seq(length=nbrOfSegments)) {
    verbose && enter(verbose, sprintf("Segment %d of %d", ss, nbrOfSegments));
    seg <- segs[ss,];

    chr <- seg[["chromosome"]];
    chrTag <- sprintf("chr%02d", chr);

    for (what in c("tcn", "dh")) {
      keys <- paste(what, c("Start", "End"), sep="");
      xStart <- seg[[keys[1]]];
      xEnd <- seg[[keys[2]]];
      verbose && printf(verbose, "[xStart,xEnd] = [%.0f,%.0f]\n", xStart, xEnd);
      # Nothing todo?
      if (is.na(xStart) && is.na(xEnd)) {
        next;
      }

      stopifnot(xStart <= xEnd);

      # (b) Identify units
      units <- which(chromosome == chr & xStart <= x & x <= xEnd);

      # (c) Adjust for missing values
      if (what == "tcn") {
        value <- CT;
      } else if (what == "dh") {
        value <- rho;
      }
      keep <- which(!is.na(value[units]));
      units <- units[keep];

      # (d) Update mean
      gamma <- mean(value[units]);

      # Sanity check
      stopifnot(length(units) == 0 || !is.na(gamma));

      # Update the segment boundaries, estimates and counts
      key <- paste(what, "mean", sep=".");
      seg[[key]] <- gamma;
    }

    segs[ss,] <- seg;

    verbose && exit(verbose);
  } # for (ss ...)

  verbose && enter(verbose, "Update (C1,C2) per segment");
  # Append (C1,C2) estimates
  tcn <- segs$tcnMean;
  dh <- segs$dhMean;
  C1 <- 1/2*(1-dh)*tcn;
  C2 <- tcn - C1;
  segs$c1Mean <- C1;
  segs$c2Mean <- C2;
  verbose && exit(verbose);

  # Return results
  res <- fit;
  res$output <- segs;

  verbose && exit(verbose);

  res;
}, private=TRUE) # updateMeans()



############################################################################
# HISTORY:
# 2011-06-14
# o Updated code to recognize new column names.
# 2011-01-18
# o BUG FIX: Fields 'tcnSegRows' and 'dhSegRows' were not updated by
#   mergeTwoSegments() for PairedPSCBS.
# 2011-01-14
# o Moved extractByRegions() and estimateStdDevForHeterozygousBAF() to
#   psCBS v0.9.36.
# o Now extractByRegions() utilizes the 'segRows' field.
# o Added estimateStdDevForHeterozygousBAF().
# 2011-01-12
# o Added updateMeans() for PairedPSCBS.
# o Added dropByRegions().
# o Added extractByRegions() and extractByRegion().
# o Created.
############################################################################
