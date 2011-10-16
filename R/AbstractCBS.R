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
# @synopsis
#
# \arguments{
#   \item{fit}{A @list structure containing the segmentation results.}
#   \item{sampleName}{A @character string.}
#   \item{...}{Not used.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
# 
# @author
#*/###########################################################################
setConstructorS3("AbstractCBS", function(fit=list(), sampleName=fit$sampleName, ...) {
  # Argument 'sampleName':
  if (!is.null(sampleName)) {
    sampleName <- Arguments$getCharacter(sampleName);
  }

  fit$sampleName <- sampleName;
  extend(fit, "AbstractCBS");
})


setMethodS3("print", "AbstractCBS", function(x, ...) {
  # To please R CMD check
  fit <- x;

  segs <- as.data.frame(fit, ...);
  print(segs);
}, protected=TRUE)




setMethodS3("all.equal", "AbstractCBS", function(target, current, check.attributes=FALSE, ...) {
  # Compare class attributes
  res <- all.equal(class(target), class(current));
  if (!isTRUE(res)) {
    return(res);
  }

  # Don't compare other attributes
  NextMethod("all.equal", target, current, check.attributes=check.attributes, ...);
}, protected=TRUE)



###########################################################################/**
# @RdocMethod getSampleName
# @aliasmethod sampleName
#
# @title "Gets the name of the sample segmented"
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
#   Returns a @character string.
# }
# 
# @author
#
# \seealso{
#   @seemethod "setSampleName".
#   @seeclass.
# }
#*/###########################################################################  
setMethodS3("getSampleName", "AbstractCBS", function(fit, ...) {
  name <- fit$sampleName;
  if (is.null(name)) {
    name <- as.character(NA);
  }
  name;
}, protected=TRUE)

setMethodS3("sampleName", "AbstractCBS", function(fit, ...) {
  getSampleName(fit);
}, protected=TRUE)


###########################################################################/**
# @RdocMethod setSampleName
# @aliasmethod sampleName<-
#
# @title "Sets the name of the sample segmented"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{name}{A @character string.}
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns (invisibly) an updated object.
# }
# 
# @author
#
# \seealso{
#   @seeclass.
# }
#
# @keyword internal
#*/###########################################################################  
setMethodS3("setSampleName", "AbstractCBS", function(fit, name, ...) {
  # Argument 'value':
  name <- Arguments$getCharacter(name);

  fit$sampleName <- name;

  invisible(fit);
}, protected=TRUE)


setMethodS3("sampleName<-", "AbstractCBS", function(x, value) {
  setSampleName(x, value);
}, protected=TRUE, addVarArgs=FALSE)

"sampleName<-" <- function(x, value) {
  UseMethod("sampleName<-");
}



###########################################################################/**
# @RdocMethod getLocusData
#
# @title "Gets the locus-level data"
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
#   Returns a JxL @data.frame, where J in the number of loci,
#   and L is the number of locus-specific fields.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################  
setMethodS3("getLocusData", "AbstractCBS", abstract=TRUE);

setMethodS3("setLocusData", "AbstractCBS", abstract=TRUE, protected=TRUE);


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
#  \item{splitters, ...}{Arguments passed to @seemethod "getLocusData".}
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
setMethodS3("nbrOfLoci", "AbstractCBS", function(fit, splitters=FALSE, ...) {
  data <- getLocusData(fit, splitters=splitters, ...);
  nrow(data);
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

setMethodS3("setSegments", "AbstractCBS", abstract=TRUE, protected=TRUE);



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
#  \item{splitters, ...}{Arguments passed to @seemethod "getSegments".}
# }
#
# \value{
#   Returns an @integer.
# }
#
# @author
#
# \seealso{
#   @seemethod "nbrOfChangePoints"
#   @seemethod "nbrOfChromosomes"
#   @seeclass
# }
#*/###########################################################################  
setMethodS3("nbrOfSegments", "AbstractCBS", function(this, splitters=FALSE, ...) {
  nrow(getSegments(this, splitters=splitters, ...));
})



###########################################################################/**
# @RdocMethod nbrOfChangePoints
#
# @title "Gets the number of change points"
#
# \description{
#   @get "title", which is defined as the number of segments minus
#   the number of chromosomes.
# }
# 
# @synopsis
#
# \arguments{
#  \item{splitters, ...}{Arguments passed to @seemethod "getSegments".}
# }
#
# \value{
#   Returns an @integer.
# }
#
# @author
#
# \seealso{
#   @seemethod "nbrOfSegments"
#   @seemethod "nbrOfChromosomes"
#   @seeclass
# }
#*/###########################################################################  
setMethodS3("nbrOfChangePoints", "AbstractCBS", function(fit, ...) {
  nbrOfSegments(fit, splitters=FALSE) - nbrOfChromosomes(fit);
})



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
#*/###########################################################################  
setMethodS3("as.data.frame", "AbstractCBS", function(x, ...) {
  getSegments(x, ...);
}, protected=TRUE)




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
#  \item{...}{Arguments passed to @seemethod "getSegments".}
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
  segs <- getSegments(this, ...);
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
setMethodS3("nbrOfChromosomes", "AbstractCBS", function(this, ...) {
  length(getChromosomes(this, ...));
})


setMethodS3("getSegmentSizes", "AbstractCBS", abstract=TRUE);

setMethodS3("extractCNs", "AbstractCBS", abstract=TRUE);

setMethodS3("sampleCNs", "AbstractCBS", function(fit, size=NULL, ...) {
  data <- extractCNs(fit, ...);

  if (!is.null(size)) {
    sizes <- getSegmentSizes(fit, ...);
    # Sanity check
    stopifnot(length(sizes) == nrow(data));
    idxs <- sample(nrow(data), size=size, replace=TRUE, prob=sizes);
    data <- data[idxs,,drop=FALSE];
  }

  data;
})

setMethodS3("updateMeans", "AbstractCBS", abstract=TRUE, protected=TRUE);


############################################################################
# HISTORY:
# 2011-10-16
# o Added sampleCNs() for AbstractCBS.
# o Added abstract getSegmentSizes() for AbstractCBS.
# o Added abstract extractCNs() for AbstractCBS.
# 2011-10-08
# o Added abstract updateMeans() for AbstractCBS.
# o Added all.equal() for AbstractCBS.
# o Added nbrOfChangePoints() for AbstractCBS.
# 2011-10-02
# o Created.
############################################################################
