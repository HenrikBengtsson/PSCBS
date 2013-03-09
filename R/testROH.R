###########################################################################/**
# @set "class=numeric"
# @RdocMethod testROH
#
# @title "Tests if a segment is in Run-of-Homozygosity (ROH)"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{betaN}{A @numeric @vector of J matched normal BAFs in [0,1]
#        (due to noise, values may be slightly outside as well) or @NA
#        for non-polymorphic loci.}
#   \item{muN}{An optional @numeric @vector of J genotype calls in 
#        \{0,1/2,1\} for AA, AB, and BB, respectively, 
#        and @NA for non-polymorphic loci.}
#   \item{csN}{A @numeric @vector of J genotype confidence scores.
#        If @NULL, ad hoc scores calculated from \code{betaN} are used.}
#   \item{minNbrOfSnps}{Minimum number of SNPs required to test segment.
#        If not tested, @NA is returned.}
#   \item{delta}{A @double scalar specifying the maximum (weighted) 
#        proportion of heterozygous SNPs allowed in an ROH region.}
#   \item{...}{Not used.}
#  \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns a @logical.
# }
#
# @author "PN, HB"
#
# @keyword internal
#*/###########################################################################
setMethodS3("testROH", "numeric", function(betaN, muN, csN=NULL, minNbrOfSnps=1, delta=1/12, ..., verbose=FALSE) {
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
  stopifnot(isZero(wnSnps - 1.0, eps=sqrt(.Machine$double.eps)));

  propHets <- wnHets/wnSnps;
  verbose && print(verbose, propHets);

  call <- (propHets < delta);
  verbose && print(verbose, call);

  # Record parameter settings
  attr(call, "minNbrOfSnps") <- minNbrOfSnps;
  attr(call, "delta") <- delta;

  verbose && exit(verbose);
  
  call;
}) # testROH()


##############################################################################
# HISTORY
# 2013-03-08 [HB]
# o Added Rdoc help.
# 2012-05-30 [HB]
# o Now testROH() return parameter settings as attributes.
# 2011-11-21 [HB]
# o BUG FIX: The internal sanity check on weights was slightly too
#   conservative.
# 2011-11-12 [HB]
# o Renamed argument 'tauROH' to 'delta', cf. how we do for AB and LOH.
# o Added argument 'minNbrOfSnps' to testROH().
# o Added verbose output.
# 2011-11-12 [PN]
# o Implemented a naive caller based on the weighted proportion of hets 
#   in the segment.
# 2011-11-04 [HB]
# o Added skeleton for testROH().
# o Created.
##############################################################################
