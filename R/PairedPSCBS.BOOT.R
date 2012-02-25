setMethodS3("bootstrapDHByRegion", "PairedPSCBS", function(fit, B=100, statsFcn=function(x) quantile(x, probs=c(0.025, 0.050, 0.95, 0.975)), by=c("betaTN", "betaT"), ..., force=FALSE, verbose=FALSE) {
  # WORKAROUND: If Hmisc is loaded after R.utils, it provides a buggy
  # capitalize() that overrides the one we want to use. Until PSCBS
  # gets a namespace, we do the following workaround. /HB 2011-07-14
  capitalize <- R.utils::capitalize;


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'B':
  B <- Arguments$getInteger(B, range=c(1,Inf));

  # Argument 'by':
  by <- match.arg(by);

  # Argument 'force':
  force <- Arguments$getLogical(force);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Resample DH signals and reestimate DH mean levels");

  data <- getLocusData(fit);
  segs <- getSegments(fit);

  # Find estimates to be done
  stats <- statsFcn(1);
  stopifnot(!is.null(names(stats)));
  statsNames <- sprintf("dh%s", capitalize(names(stats)));
  isDone <- is.element(statsNames, names(segs));

  # Already done?
  if (!force && all(isDone)) {
    verbose && cat(verbose, "Already done. Skipping.");
    verbose && exit(verbose);
    return(fit);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Precalculate signals
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get (x,BAF) data
  chromosome <- data$chromosome;
  x <- data$x;
  muN <- data$muN;
  beta <- data[[by]];

  # SNPs are identifies as those loci that have non-missing 'betaTN' & 'muN'
  isSnp <- (!is.na(beta) & !is.na(muN));
  nbrOfSnps <- sum(isSnp);
  verbose && cat(verbose, "Number of SNPs: ", nbrOfSnps);

  # DH is by definition only defined for heterozygous SNPs.  
  # For simplicity, we set it to be NA for non-heterozygous loci.
  isHet <- isSnp & (muN == 1/2);
  verbose && printf(verbose, "Number of heterozygous SNPs: %d (%.2f%%)\n",
                                     sum(isHet), 100*sum(isHet)/nbrOfSnps);

  # Extract subset for heterozygous SNPs
  chromosome <- chromosome[isHet];
  x <- x[isHet];
  beta <- beta[isHet];

  # Calculate DHs for heterozygous SNPs
  rho <- 2*abs(beta - 1/2);

  # Drop missing values  (Is that ok for bootstrapping? /HB 2010-09-16)
  keep <- (is.finite(chromosome) & is.finite(x));
  chromosome <- chromosome[keep];
  x <- x[keep];
  rho <- rho[keep];

  # Not needed anymore
  rm(beta, muN, isHet, keep);

  # Sanity checks
  stopifnot(all(is.finite(chromosome)));
  stopifnot(all(is.finite(x)));
  stopifnot(all(is.finite(rho)));

  listOfDhLociNotPartOfSegment <- fit$listOfDhLociNotPartOfSegment;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Resample DH within segments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  nbrOfSegments <- nrow(segs);
  naValue <- as.double(NA);

  rhoMeanList <- vector("list", nbrOfSegments);
  for (jj in seq(length=nbrOfSegments)) {
    verbose && enter(verbose, sprintf("DH segment #%d of %d", jj, nbrOfSegments));
    segJJ <- segs[jj,,drop=FALSE];

    # Identify loci in segment
    chr <- segJJ$chromosome[1];
    start <- segJJ$dhStart[1];
    stop <- segJJ$dhEnd[1];
    nbrOfDHs <- segJJ[,"dhNbrOfLoci"];
    if (is.na(nbrOfDHs)) nbrOfDHs <- 0L;

    units <- whichVector(chr == chromosome & start <= x & x <= stop);
    nbrOfUnits <- length(units);

    # Special case?
    if (nbrOfUnits > nbrOfDHs) {
      verbose && cat(verbose, "All loci in DH segment: ", nbrOfUnits);
      verbose && cat(verbose, "Used loci in DH segment: ", nbrOfDHs);

      # Sanity check
      stopifnot(!is.null(listOfDhLociNotPartOfSegment));

      tcnId <- segJJ[,"tcnId"];
      dhId <- segJJ[,"dhId"];
      dhLociNotPartOfSegment <- listOfDhLociNotPartOfSegment[[tcnId]];
      # Sanity check
      stopifnot(!is.null(dhLociNotPartOfSegment));

      lociToExclude <- dhLociNotPartOfSegment[[dhId]];
      verbose && cat(verbose, "Excluding loci that belongs to a flanking segment: ", length(lociToExclude));
      units <- setdiff(units, lociToExclude);
      nbrOfUnits <- length(units);
    }

    # Sanity check
    stopifnot(nbrOfUnits == nbrOfDHs);

    if (nbrOfUnits >= 1) {
      # Sanity check
      mu <- mean(rho[units], na.rm=FALSE);
      dMu <- (mu - segJJ$dhMean);
      tol <- 0.0005;
      if (abs(dMu) > tol) {
        str(list(nbrOfUnits=nbrOfUnits, dhNbrOfLoci=segJJ$dhNbrOfLoci, mu=mu, dhMean=segJJ$dhMean, dMu=dMu, "abs(dMu)"=abs(dMu), "min(x[units])"=min(x[units])));
#        stop("INTERNAL ERROR");
      }

      # Bootstrap B times
      rhoMean <- rep(naValue, times=B);
      for (bb in seq(length=B)) {
        # Resample indices
        unitsS <- resample(units, size=nbrOfUnits, replace=TRUE);
  
        # Resample data
        rhoB <- rho[unitsS];

        # Calculate new mean level (no NAs here)
        rhoMean[bb] <- mean(rhoB, na.rm=FALSE);
      } # for (bb ...)
    } else {
      rhoMean <- double(0);
    }
    rhoMeanList[[jj]] <- rhoMean;

    verbose && exit(verbose);
  } # for (jj ...)

  # Not needed anymore
  rm(x, rho);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Calculate bootstrap statistics
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  verbose && enter(verbose, "Calculating bootstrap statistics");
  statsList <- vector("list", length(rhoMeanList));
  for (jj in seq(along=rhoMeanList)) {
    verbose && enter(verbose, sprintf("DH segment #%d of %d", 
                                             jj, length(rhoMeanList)));
    # Get bootstrap sample
    rhoMean <- rhoMeanList[[jj]];

    # Calculate statistics
    stats <- statsFcn(rhoMean);
    # Store
    statsList[[jj]] <- stats;

    verbose && exit(verbose);
  } # for (jj ...)
  verbose && exit(verbose);

  # Flatten statistics, if possible
  stats <- sapply(statsList, FUN=function(x) x);
  if (!is.matrix(stats)) {
    stats <- as.matrix(stats);
  } else {
    stats <- t(stats);
  }

  # Sanity check
  stopifnot(is.matrix(stats));

  # Column matrix
  colnames(stats) <- statsNames;

  # Store results
  segs <- cbind(segs, stats);

  # Drop previously estimated values
  dups <- duplicated(colnames(segs), fromLast=TRUE);
  if (any(dups)) {
    stats <- stats[,!dups, drop=FALSE];
  }
  
  fitB <- fit;
  fitB$output <- segs;

  verbose && exit(verbose);

  fitB;
}, deprecated=TRUE) # bootstrapDHByRegion()



##############################################################################
# HISTORY
# 2012-02-24
# o Added argument 'force' to bootstrapDHByRegion().
# 2011-06-14
# o Updated code to recognize new column names.
# 2010-11-22
# o DEPRECATED: bootstrapDHByRegion() should no longer be used.
# 2010-11-03 [HB]
# o ROBUSTNESS: Now bootstrapDHByRegion() uses resample() of R.utils.
# 2010-11-01 [HB]
# o Now bootstrapDHByRegion() estimates more quantiles.
# o BUG FIX: bootstrapDHByRegion() would give an error if only a single
#   quantile was requested.
# o BUG FIX: bootstrapDHByRegion() would give "Error in if (nbrOfUnits > 
#   segJJ[, "dh.num.mark"]) { : missing value where TRUE/FALSE needed" when
#   'dh.num.mark' was NA.
# 2010-10-25 [HB]
# o BUG FIX: Now bootstrapDHByRegion() for PairedPSCBS handles the rare case
#   when markers with the same positions are split in two different segments.
# o BUG FIX: bootstrapDHByRegion() for PairedPSCBS would bootstrap from the
#   incorrect set of loci when the DH region contained only one locus.
# o BUG FIX: bootstrapDHByRegion() for PairedPSCBS would bootstrap from the
#   incorrect set of loci if more than one chromosome was available.
# 2010-09-16 [HB]
# o Added bootstrapDHByRegion(), which is what is used by paired PSCBS.
# o Created.
##############################################################################
