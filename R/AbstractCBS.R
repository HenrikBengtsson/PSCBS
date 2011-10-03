###########################################################################/**
# @RdocClass AbstractCBS
#
# @title "The AbstractCBS class"
#
# \description{
#  @classhierarchy
#
#  All CBS-style segmentation results extend this class, e.g.
#  @see "CBS" and @see "PairedPSCBS".
# }
# 
# \usage{AbstractCBS(fit=list(), ...)}
#
# \arguments{
#   \item{fit}{A @list structure containing the segmentation results.}
#   \item{...}{Not used.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
# 
# @author
#*/###########################################################################
setConstructorS3("AbstractCBS", function(fit=list(), ...) {
  extend(fit, "AbstractCBS");
})


setMethodS3("print", "AbstractCBS", function(x, ...) {
  # To please R CMD check
  fit <- x;

  segs <- as.data.frame(fit, ...);
  print(segs);
}, private=TRUE)


setMethodS3("getSampleName", "AbstractCBS", function(fit, ...) {
  name <- fit$sampleName;
  if (is.null(name)) {
    name <- as.character(NA);
  }
  name;
}, private=TRUE)

setMethodS3("sampleName", "AbstractCBS", function(fit, ...) {
  getSampleName(fit);
}, private=TRUE)

"sampleName<-" <- function(x, value) {
  UseMethod("sampleName<-");
}

setMethodS3("sampleName<-", "AbstractCBS", function(x, value) {
  fit <- x;

  # Argument 'value':
  value <- Arguments$getCharacter(value);

  fit$sampleName <- value;
  fit;
}, private=TRUE, addVarArgs=FALSE)


###########################################################################/**
# @RdocMethod as.data.frame
#
# @title "Gets the table of segments"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns a @data.frame, where each row corresponds to 
#   a unique segment.
# }
# 
# @author
#
# \seealso{
#   Utilizes @seemethod "getSegments".
#   @seeclass.
# }
#
# @keyword internal
#*/###########################################################################  
setMethodS3("as.data.frame", "AbstractCBS", function(x, ...) {
  getSegments(x, ...);
})




###########################################################################/**
# @RdocMethod getChromosomes
#
# @title "Gets the set of chromosomes"
#
# \description{
#   @get "title" in the segmentation result.
# }
# 
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a unique and sorted @vector of chromosomes segmented.
# }
#
# @author
#
# \seealso{
#   @seemethod "nbrOfChromosomes".
#   @seeclass
# }
#*/###########################################################################  
setMethodS3("getChromosomes", "AbstractCBS", function(this, ...) {
  segs <- getSegments(this);
  chromosomes <- sort(unique(segs$chromosome), na.last=TRUE);

  # Drop NA dividers
  if (length(chromosomes) > 1) {
    chromosomes <- chromosomes[!is.na(chromosomes)];
  }

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
#   @seemethod "getChromosomes".
#   @seeclass
# }
#*/###########################################################################  
setMethodS3("nbrOfChromosomes", "AbstractCBS", function(this, ...) {
  length(getChromosomes(this, ...));
})



###########################################################################/**
# @RdocMethod getSegments
#
# @title "Gets the segments"
#
# \description{
#   @get "title".
# }
# 
# @synopsis
#
# \arguments{
#  \item{splitters}{If @TRUE, "splitters" between chromosomes are 
#     preserved, otherwise dropped.}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a SxK @data.frame, where S in the number of segments,
#   and K is the number of segment-specific fields.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################  
setMethodS3("getSegments", "AbstractCBS", abstract=TRUE);



###########################################################################/**
# @RdocMethod nbrOfSegments
#
# @title "Gets the number of segments"
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
setMethodS3("nbrOfSegments", "AbstractCBS", function(this, splitters=FALSE, ...) {
  nrow(getSegments(this, splitters=splitters, ...));
})



setMethodS3("extractByChromosome", "AbstractCBS", function(x, chromosome, ...) {
  # To please R CMD check
  this <- x;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'chromosome':
  chromosome <- Arguments$getInteger(chromosome, disallow=c("NaN", "Inf"));

  extractByChromosomes(this, chromosomes=chromosome, ...);
}, protected=TRUE)



setMethodS3("extractByChromosomes", "AbstractCBS", abstract=TRUE, protected=TRUE);




###########################################################################/**
# @set "class=PSCBS"
# @RdocMethod append
#
# @title "Appends one segmentation result to another"
#
# \description{
#   @get "title".
# }
# 
# @synopsis
#
# \arguments{
#  \item{x, other}{The two @see "AbstractCBS" objects to be combined.}
#  \item{addSplit}{If @TRUE, a "divider" is added between chromosomes.}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a object of the same class as argument \code{x}.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################  
setMethodS3("append", "AbstractCBS", abstract=TRUE);


setMethodS3("plotTracks", "AbstractCBS", abstract=TRUE);


############################################################################
# HISTORY:
# 2011-10-02
# o Created.
############################################################################
