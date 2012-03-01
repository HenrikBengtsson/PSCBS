###########################################################################/**
# @RdocClass AbstractCNData
#
# @title "The AbstractCNData class"
#
# \description{
#  @classhierarchy
#
#  An AbstractCNData object holds copy number data.
# }
# 
# @synopsis
#
# \arguments{
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
setConstructorS3("AbstractCNData", function(chromosome=NULL, x=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!is.null(chromosome)) {
    # Argument 'chromosome':
    if (is.data.frame(chromosome)) {
      data <- chromosome;
      chromosome <- data$chromosome;
      x <- data$x;
    }
  
    # Argument 'chromosome':
    disallow <- c("Inf");
    chromosome <- Arguments$getIntegers(chromosome, range=c(0,Inf), disallow=disallow);
    nbrOfLoci <- length(chromosome);
    length2 <- rep(nbrOfLoci, times=2);
  
    # Argument 'x':
    if (is.null(x)) {
      x <- seq(length=nbrOfLoci);
    } else {
      disallow <- c("Inf");
      x <- Arguments$getDoubles(x, length=length2, disallow=disallow);
    }
   
    # Build data.frame()
    data <- data.frame(chromosome=chromosome, x=x);
  
    # Any extra locus signals?
    args <- list(...);
    keep <- !sapply(args, FUN=is.null);
    args <- args[keep];
    if (length(args) > 0) {
      data <- cbind(data, args);
    }
  } else {
    # Dummy
    data <- data.frame();
  }

  extend(data, "AbstractCNData");
})


setMethodS3("getLocusData", "AbstractCNData", function(this, ...) {
  data <- this;
  class(data) <- c("data.frame");
  data;
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
setMethodS3("nbrOfLoci", "AbstractCNData", function(this, ...) {
  data <- getLocusData(this, ...);
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
setMethodS3("getChromosomes", "AbstractCNData", function(this, ...) {
  data <- getLocusData(this, ...);
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
setMethodS3("nbrOfChromosomes", "AbstractCNData", function(this, ...) {
  length(getChromosomes(this, ...));
})


setMethodS3("hasKnownPositions", "AbstractCNData", function(this, ...) {
  data <- getLocusData(this, ...);
  ok <- (!is.na(data$chromosome) & !is.na(data$x));
  ok;
}, protected=TRUE)



setMethodS3("orderAlongGenome", "AbstractCNData", function(this, ...) {
  data <- getLocusData(this, ...);

  o <- order(data$chromosome, data$x, decreasing=FALSE, na.last=TRUE);
  rm(data); # Not needed anymore

  # Any change?
  if (any(o != seq(along=o))) {
    this <- this[o,,drop=FALSE];
  }
  rm(o); # Not needed anymore

  this;
}, protected=TRUE)



setMethodS3("findLargeGaps", "AbstractCNData", function(chromosome, ...) {
  # To please R CMD check
  this <- chromosome;

  data <- getLocusData(this, ...);
  findLargeGaps(chromosome=data$chromosome, x=data$x, ...); 
}) # findLargeGaps()



############################################################################
# HISTORY:
# 2012-02-29
# o Added getLocusData() for AbstractCNData.
# o Added findLargeGaps() for AbstractCNData.
# o Created.
############################################################################
