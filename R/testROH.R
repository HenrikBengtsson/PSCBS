#
# \arguments{
#   \item{tauROH}{A @double in [0,1] specifying the maximum fraction of hets
#     for the segment to be called "ROH". Defaults to 1/12 (which is much
#     smaller than the typical fraction of hets.}
# }
#
# \value{
#  Returns @TRUE if segment is ROH, otherwise @FALSE.
#  If it was not possible to make a call, then @NA is returned.
# }
setMethodS3("testROH", "numeric", function(betaN, muN, csN=NULL, minNbrOfSnps=1, tauROH=1/12, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'betaN':
  betaN <- Arguments$getDoubles(betaN);
  nbrOfSNPs <- length(betaN);

  # Argument 'muN':
  length2 <- rep(nbrOfSNPs, times=2);
  muN <- Arguments$getDoubles(muN, range=c(0,1), length=length2);

  # Argument 'csN':
  if (!is.null(csN)) {
    csN <- Arguments$getDoubles(csN, range=c(0,1), length=length2);
  }

  # Argument 'minNbrOfSnps':
  minNbrOfSnps <- Arguments$getInteger(minNbrOfSnps, range=c(1,Inf));

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Testing for ROH");

  # Default ROH call
  call <- NA;

  nbrOfSnps <- length(betaN);
  verbose && cat(verbose, "Number of SNPs: ", nbrOfSnps);

  # Nothing todo?
  if (nbrOfSnps < minNbrOfSnps) {
    verbose && exit(verbose);
    return(call);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Calculate genotype confidence scores?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (is.null(csN)) {
    verbose && enter(verbose, "Calculating confidence scores");
    # Assuming naive genotyping a'la aroma.light::callNaiveGenotypes()
    # was used to call genotypes 'muN' from 'betaN'.
    
    # AD HOC: We also have to assume that the thresholds were 1/3 and 2/3.
    a <- 1/3;  # was fit$x[1];
    b <- 2/3;  # was fit$x[2];

    # AD HOC: We have to make some assumption about which SNPs are diploid.
    # Assume all for now
    isDiploid <- rep(TRUE, times=nbrOfSNPs);

    # KNOWN ISSUE: Scores for homozygotes are in [0,1/3], whereas
    # heterzygotes are in [0,1/6]. /PN 2011-11-11
    csN[isDiploid] <- rowMins(abs(cbind(betaN[isDiploid]-a, betaN[isDiploid]-b)));
    verbose && exit(verbose);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Call ROH
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # 0-1 weights (just to make sure)
  # Weights summing to one
  w <- csN / sum(csN, na.rm=TRUE);

  # Identify heterozygous SNPs
  isHet <- (muN == 1/2);
  verbose && print(verbose, summary(isHet));

  wnHets <- sum(isHet*w, na.rm=TRUE);
  wnSnps <- sum(w, na.rm=TRUE);  # == 1 /HB

  # Sanity check
  stopifnot(isZero(wnSnps - 1.0));

  propHets <- wnHets/wnSnps;
  verbose && print(verbose, propHets);

  call <- (propHets < tauROH);
  verbose && print(verbose, call);

  verbose && exit(verbose);
  
  call;
}) # testROH()


##############################################################################
# HISTORY
# 2011-11-12 [HB]
# o Added argument 'minNbrOfSnps' to testROH().
# o Added verbose output.
# 2011-11-12 [PN]
# o Implemented a naive caller based on the weighted proportion of hets 
#   in the segment.
# 2011-11-04 [HB]
# o Added skeleton for testROH().
# o Created.
##############################################################################
