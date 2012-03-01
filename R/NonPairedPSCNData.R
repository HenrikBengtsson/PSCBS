###########################################################################/**
# @RdocClass NonPairedPSCNData
#
# @title "The NonPairedPSCNData class"
#
# \description{
#  @classhierarchy
#
#  A NonPairedPSCNData object holds parent-specific copy number data.
#  Two NonPairedPSCNData objects for a matched tumor-normal pair can
#  be combined into a @see "PairedPSCNData" object.
# }
# 
# @synopsis
#
# \arguments{
#   \item{C}{A @numeric @vector of J tumor total copy number (TCN)
#        ratios in [0,+@Inf) (due to noise, small negative values are
#        also allowed).  The TCN ratios are typically scaled such that
#        copy-neutral diploid loci have a mean of two.}
#   \item{beta}{A @numeric @vector of J tumor allele B fractions (BAFs)
#        in [0,1] (due to noise, values may be slightly outside as well)
#        or @NA for non-polymorphic loci.}
#   \item{mu}{An optional @numeric @vector of J genotype calls in 
#        \{0,1/2,1\} for AA, AB, and BB, respectively, 
#        and @NA for non-polymorphic loci.
#        If not given, they are estimated from the normal BAFs using
#        @see "aroma.light::callNaiveGenotypes" as described in [2].}
#   \item{isSNP}{An optional @logical @vector of length J specifying
#        whether each locus is a SNP or not (non-polymorphic loci).}
#   \item{chromosome}{(Optional) An @integer scalar (or a @vector of length J),
#        which can be used to specify which chromosome each locus belongs to
#        in case multiple chromosomes are segments.
#        This argument is also used for annotation purposes.}
#   \item{x}{Optional @numeric @vector of J genomic locations.
#            If @NULL, index locations \code{1:J} are used.}
#   \item{...}{Optional named locus-specific signal @vectors of length J.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
#*/###########################################################################
setConstructorS3("NonPairedPSCNData", function(chromosome=NULL, x=NULL, isSNP=NULL, mu=NULL, C=NULL, beta=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!is.null(chromosome)) {
    # Is first argument special?
    if (is.data.frame(chromosome)) {
      data <- chromosome;

      chromosome <- data$chromosome;
      x <- data$x;
      isSNP <- data$isSNP;
      mu <- data$mu;
      C <- data$C; 
      beta <- data$beta;
    }

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Required arguments
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Argument 'C':
    disallow <- c("Inf");
    C <- Arguments$getDoubles(C, disallow=disallow);
  
    nbrOfLoci <- length(C);
    length2 <- rep(nbrOfLoci, times=2);
  

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Mutually optional arguments (that are not validated in superclass)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Argument 'beta':
    if (!is.null(beta)) {
      beta <- Arguments$getDoubles(beta, length=length2, disallow="Inf");
    }
  
    if (is.null(beta) && is.null(mu)) {
      throw("If argument 'beta' is not given, then 'mu' must be.");
    }
  

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Optional arguments (that are not validated in superclass)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Argument 'chromosome':
    if (is.null(chromosome)) {
      chromosome <- 0L; 
    } else {
      disallow <- c("Inf");
      chromosome <- Arguments$getIntegers(chromosome, range=c(0,Inf), disallow=disallow);
    }
  
    if (length(chromosome) == 1) {
      chromosome <- rep(chromosome, times=nbrOfLoci);
    } else {
      chromosome <- Arguments$getIntegers(chromosome, length=length2, disallow=disallow);
    }
  }

  extend(AbstractPSCNData(chromosome=chromosome, x=x, isSNP=isSNP, mu=mu, C=C, beta=beta, ...), "NonPairedPSCNData");
})


setMethodS3("as.NonPairedPSCNData", "NonPairedPSCNData", function(this, ...) {
  this;
})

setMethodS3("as.NonPairedPSCNData", "data.frame", function(this, ...) {
  data <- this;
  NonPairedPSCNData(data, ...);
})



setMethodS3("callNaiveGenotypes", "NonPairedPSCNData", function(this, censorAt=c(0,1), force=FALSE, ..., verbose=FALSE) {
  data <- this;

  require("aroma.light") || throw("Package not loaded: aroma.light");

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Calling genotypes from allele B fractions (BAFs)");
  mu <- data$mu;
  if (!force && !is.null(mu)) {
    verbose && cat(verbose, "Already done. Skipping.");
    verbose && exit(verbose);
    return(data);
  }

  # Call genotypes from BAFs
  verbose && cat(verbose, "Allele B fractions (BAFs):");
  verbose && str(verbose, data$beta);

  mu <- callNaiveGenotypes(data$beta, censorAt=censorAt, ...);
  verbose && cat(verbose, "Called genotypes:");
  verbose && str(verbose, mu);
  verbose && print(verbose, table(mu));

  data$mu <- mu;
  verbose && exit(verbose);

  data;
})



setMethodS3("callSegmentationOutliers", "NonPairedPSCNData", function(y, ...) {
  # To please R CMD check
  data <- y;

  C <- data$C;
  callSegmentationOutliers(y=C, chromosome=data$chromosome, x=data$x, ...);
}) # callSegmentationOutliers()


setMethodS3("dropSegmentationOutliers", "NonPairedPSCNData", function(C, ...) {
  # To please R CMD check
  data <- C;

  isOutlier <- callSegmentationOutliers(data, ...);
  naValue <- as.double(NA);
  data$C[isOutlier] <- naValue;
  data;
}) # dropSegmentationOutliers()



############################################################################
# HISTORY:
# 2012-02-29
# o Added callSegmentationOutliers() and dropSegmentationOutliers() 
#   for NonPairedPSCNData.
# o Added as.NonPairedPSCNData().
# o Created.
############################################################################
