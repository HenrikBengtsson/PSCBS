###########################################################################/**
# @RdocClass PairedPSCNData
#
# @title "The PairedPSCNData class"
#
# \description{
#  @classhierarchy
#
#  A PairedPSCNData object holds paired tumor-normal parent-specific 
#  copy number data.
# }
# 
# @synopsis
#
# \arguments{
#   \item{CT}{A @numeric @vector of J tumor total tumor copy number (TCN)
#        ratios in [0,+@Inf) (due to noise, small negative values are
#        also allowed).  The TCN ratios are typically scaled such that
#        copy-neutral diploid loci have a mean of two.}
#   \item{betaT}{A @numeric @vector of J tumor allele B fractions (BAFs)
#        in [0,1] (due to noise, values may be slightly outside as well)
#        or @NA for non-polymorphic loci.}
#   \item{betaN}{A @numeric @vector of J matched normal BAFs in [0,1]
#        (due to noise, values may be slightly outside as well) or @NA
#        for non-polymorphic loci.}
#   \item{muN}{An optional @numeric @vector of J genotype calls in 
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
setConstructorS3("PairedPSCNData", function(CT=NULL, betaT, betaN=NULL, muN=NULL, isSNP=NULL, chromosome=0, x=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!is.null(CT)) {
    # Argument 'CT':
    if (is.data.frame(CT)) {
      data <- CT;
      CT <- data$CT; 
      betaT <- data$betaT;
      betaN <- data$betaN;
      muN <- data$muN;
      isSNP <- data$isSNP;
      chromosome <- data$chromosome;
      x <- data$x;
    }
    disallow <- c("Inf");
    CT <- Arguments$getDoubles(CT, disallow=disallow);
  
    nbrOfLoci <- length(CT);
    length2 <- rep(nbrOfLoci, times=2);
  
    # Argument 'betaT':
    betaT <- Arguments$getDoubles(betaT, length=length2, disallow="Inf");
  
    # Argument 'betaN':
    if (!is.null(betaN)) {
      betaN <- Arguments$getDoubles(betaN, length=length2, disallow="Inf");
    }
  
    # Argument 'muN':
    if (!is.null(muN)) {
      muN <- Arguments$getDoubles(muN, length=length2, range=c(0,1), disallow="Inf");
    }
  
    if (is.null(betaN) && is.null(muN)) {
      throw("If argument 'betaN' is not given, then 'muN' must be.");
    }
  
    # Argument 'isSNP':
    if (!is.null(isSNP)) {
      isSNP <- Arguments$getLogicals(isSNP, length=length2, disallow="NA");
    }
  
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
  
    # Argument 'x':
    if (is.null(x)) {
      x <- seq(length=nbrOfLoci);
    } else {
      disallow <- c("Inf");
      x <- Arguments$getDoubles(x, length=length2, disallow=disallow);
    }
   
  
    # Build data.frame()
    data <- data.frame(chromosome=chromosome, x=x, CT=CT, betaT=betaT);
  
    # Add the optional locus signals
    if (!is.null(betaN)) {
      data$betaN <- betaN;
    }
    if (!is.null(muN)) {
      data$muN <- muN;
    }
    if (!is.null(isSNP)) {
      data$isSNP <- isSNP;
    }
  
    # Any extra locus signals?
    extra <- data.frame(...);
    if (ncol(extra) > 0) {
      data <- cbind(data, extra);
    }
  } else {
    # Dummy
    data <- data.frame();
  }

  extend(data, "PairedPSCNData");
})


setMethodS3("as.PairedPSCNData", "PairedPSCNData", function(this, ...) {
  this;
})

setMethodS3("as.PairedPSCNData", "data.frame", function(this, ...) {
  data <- this;
  PairedPSCNData(data, ...);
})



###########################################################################/**
# @RdocMethod nbrOfLoci
#
# @title "Gets the number of loci"
#
# \description{
#   @get "title".
# }
# 
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns an @integer.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################  
setMethodS3("nbrOfLoci", "PairedPSCNData", function(this, ...) {
  data <- this;
  nrow(data);
})




###########################################################################/**
# @RdocMethod getChromosomes
#
# @title "Gets the set of chromosomes"
#
# \description{
#   @get "title".
# }
# 
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a unique and sorted @vector of chromosome indices.
# }
#
# @author
#
# \seealso{
#   @seemethod "nbrOfChromosomes".
#   @seeclass
# }
#*/###########################################################################  
setMethodS3("getChromosomes", "PairedPSCNData", function(this, ...) {
  data <- this;
  chromosomes <- sort(unique(data$chromosome), na.last=TRUE);
  chromosomes;
})


###########################################################################/**
# @RdocMethod nbrOfChromosomes
#
# @title "Gets the number of chromosomes"
#
# \description{
#   @get "title".
# }
# 
# @synopsis
#
# \arguments{
#  \item{...}{Arguments passed to @seemethod "getChromosomes".}
# }
#
# \value{
#   Returns an @integer.
# }
#
# @author
#
# \seealso{
#   @seemethod "getChromosomes".
#   @seeclass
# }
#*/###########################################################################  
setMethodS3("nbrOfChromosomes", "PairedPSCNData", function(this, ...) {
  length(getChromosomes(this, ...));
})



setMethodS3("callNaiveGenotypes", "PairedPSCNData", function(this, censorAt=c(0,1), force=FALSE, ..., verbose=FALSE) {
  data <- this;

  require("aroma.light") || throw("Package not loaded: aroma.light");

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Calling genotypes from normal allele B fractions");

  muN <- data$muN;
  if (!force && !is.null(muN)) {
    verbose && cat(verbose, "Already done. Skipping.");
    verbose && exit(verbose);
    return(data);
  }

  # If muN is missing, call genotypes from betaN
  verbose && cat(verbose, "Normal (germline) BAFs:");
  verbose && str(verbose, data$betaN);

  muN <- callNaiveGenotypes(data$betaN, censorAt=censorAt, ...);
  verbose && cat(verbose, "Called genotypes:");
  verbose && str(verbose, muN);
  verbose && print(verbose, table(muN));

  data$muN <- muN;
  verbose && exit(verbose);

  data;
})


setMethodS3("normalizeTumorBoost", "PairedPSCNData", function(this, preserveScale=TRUE, force=FALSE, ..., verbose=FALSE) {
  data <- this;

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Normalizing tumor BAFs using normal BAFs (TumorBoost)");

  betaTN <- data$betaTN;
  if (!force && !is.null(betaTN)) {
    verbose && cat(verbose, "Already done. Skipping.");
    verbose && exit(verbose);
    return(data);
  }

  betaT <- data$betaT;
  betaN <- data$betaN;
  muN <- data$muN;

  betaTN <- aroma.light::normalizeTumorBoost(betaT=betaT, betaN=betaN, muN=muN, preserveScale=preserveScale);
  verbose && cat(verbose, "Normalized tumor BAFs:");
  verbose && str(verbose, betaTN);

  # Assert that no missing values where introduced
  keep <- (is.finite(betaT) & is.finite(betaN) & is.finite(muN));
  if (any(is.na(betaTN[keep]))) {
    throw("Internal error: normalizeTumorBoost() introduced missing values.");
  }
  rm(keep);

  data$betaTN <- betaTN;
  verbose && exit(verbose);

  data;
}) # normalizeTumorBoost()


setMethodS3("isSNP", "PairedPSCNData", function(this, method=c("some", "all"), force=FALSE, ...) {
  data <- this;

  # Argument 'method':
  method <- match.arg(method);

  isSNP <- data$isSNP;
  if (!force && !is.null(isSNP)) {
    verbose && cat(verbose, "Already done. Skipping.");
    verbose && exit(verbose);
    return(isSNP);
  }


  verbose && enter(verbose, "Inferring SNPs from BAFs and genotype calls");
  nbrOfLoci <- nrow(data);

  # First, assume none of the loci are SNPs
  if (method == "some") {
    isSNP <- rep(FALSE, times=nbrOfLoci);
    fcn <- function(x,y) { x | y };
  } else if (method == "all") {
    isSNP <- rep(TRUE, times=nbrOfLoci);
    fcn <- function(x,y) { x & y };
  }

  # Then remove those for which we would expect to have SNP signals
  fields <- c("muN", "betaN", "betaT", "betaTN");
  for (ff in seq(along=fields)) {
    field <- fields[ff];
    signals <- data[[field]];
    if (is.null(signals)) {
      next;
    }
    isSNP <- fcn(isSNP, !is.na(signals));
  } # for (ff ...)

  # Sanity check
  stopifnot(all(is.finite(isSNP)));

  nbrOfSNPs <- sum(isSNP);
  verbose && printf(verbose, "Number of SNPs: %d (%.4g%% of %d loci)\n", nbrOfSNPs, 100*nbrOfSNPs/nbrOfLoci, nbrOfLoci); 

  data$isSNP <- isSNP;
  verbose && exit(verbose);

  data;
}, protected=TRUE)


setMethodS3("hasKnownPositions", "PairedPSCNData", function(this, ...) {
  # To please R CMD check
  data <- this;

  ok <- (!is.na(data$chromosome) & !is.na(data$x));
  ok;
}, protected=TRUE)



setMethodS3("orderAlongGenome", "PairedPSCNData", function(this, ...) {
  # To please R CMD check
  data <- this;

  o <- order(data$chromosome, data$x, decreasing=FALSE, na.last=TRUE);
  # Any change?
  if (any(o != seq(along=o))) {
    data <- data[o,,drop=FALSE];
  }
  rm(o); # Not needed anymore
  data;
}, protected=TRUE)


setMethodS3("callSegmentationOutliers", "PairedPSCNData", function(y, ...) {
  # To please R CMD check
  data <- y;

  CT <- data$CT;
  callSegmentationOutliers(y=CT, chromosome=data$chromosome, x=data$x, ...);
}) # callSegmentationOutliers()


setMethodS3("dropSegmentationOutliers", "PairedPSCNData", function(CT, ...) {
  # To please R CMD check
  data <- CT;

  isOutlier <- callSegmentationOutliers(data, ...);
  naValue <- as.double(NA);
  data$CT[isOutlier] <- naValue;
  data;
}) # dropSegmentationOutliers()


setMethodS3("findLargeGaps", "PairedPSCNData", function(chromosome, ...) {
  # To please R CMD check
  data <- chromosome;

  findLargeGaps(chromosome=data$chromosome, x=data$x, ...); 
}) # findLargeGaps()


setMethodS3("segmentByPairedPSCBS", "PairedPSCNData", function(CT, ...) {
  # To please R CMD check
  data <- CT;

  segmentByPairedPSCBS(CT=data$CT, betaT=data$betaT, betaN=data$betaN, 
                   muN=data$muN, chromosome=data$chromosome, x=data$x, ...); 
}) # segmentByPairedPSCBS()



############################################################################
# HISTORY:
# 2012-02-29
# o Added findLargeGaps() for PairedPSCNData.
# o Added callSegmentationOutliers() and dropSegmentationOutliers() 
#   for PairedPSCNData.
# o Added segmentByPairedPSCBS() for PairedPSCNData.
# o Added as.PairedPSCNData().
# o Added isSNP() for PairedPSCNData.
# o Added callNaiveGenotypes() for PairedPSCNData.
# o Added normalizeTumorBoost() for PairedPSCNData.
# o Created.
############################################################################
