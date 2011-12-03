###########################################################################/**
# @RdocClass PSCBS
#
# @title "The PSCBS class"
#
# \description{
#  @classhierarchy
#
#  A PSCBS is an object containing results from parent-specific copy-number
#  (PSCN) segmentation.
# }
# 
# \usage{PSCBS(fit=list(), ...)}
#
# \arguments{
#   \item{fit}{A @list structure containing the PSCN segmentation results.}
#   \item{...}{Not used.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
# 
# @author
#
# \seealso{
#   @see "PairedPSCBS".
# }
#*/###########################################################################
setConstructorS3("PSCBS", function(fit=list(), ...) {
  # Argument 'fit':
  if (!is.list(fit)) {
    throw("Argument 'fit' is not a list: ", class(fit)[1]);
  }

  extend(AbstractCBS(fit, ...), "PSCBS");
})


setMethodS3("as.data.frame", "PSCBS", function(x, ...) {
  getSegments(x, splitter=TRUE, ...);
}, protected=TRUE)


setMethodS3("getLocusData", "PSCBS", function(fit, ...) {
  data <- fit$data;
  data;
}, protected=TRUE)


setMethodS3("isSegmentSplitter", "PSCBS", function(fit, ...) {
  segs <- fit$output;

  isSplitter <- lapply(segs[-1], FUN=is.na);
  isSplitter <- Reduce("&", isSplitter);

  isSplitter;
}, protected=TRUE)

 
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
#  \item{simplify}{If @TRUE, redundant and intermediate information is dropped.}#  \item{splitters}{If @TRUE, "splitters" between chromosomes are 
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
setMethodS3("getSegments", "PSCBS", function(fit, simplify=FALSE, splitters=TRUE, ...) {
  # Argument 'splitters':
  splitters <- Arguments$getLogical(splitters);

  segs <- fit$output;

  # Drop chromosome splitters?
  if (!splitters) {
    isSplitter <- isSegmentSplitter(fit);
    segs <- segs[!isSplitter,];
  }

##  if (nrow(segs) > 0) {
##    segs$id <- getSampleName(fit);
##  }

  if (simplify) {
    # If joinSegments was used (i.e. (start,end) are equal for TCN and DH)...
    if (fit$params$joinSegments) {
      # Sanity check
      stopifnot(all(segs$tcnStart == segs$dhStart, na.rm=TRUE));
      stopifnot(all(segs$tcnEnd == segs$dhEnd, na.rm=TRUE));

      names <- colnames(segs);
      keep <- !is.element(names, c("dhStart", "dhEnd"));
      segs <- segs[,keep];
      names <- colnames(segs);
      names[names == "tcnStart"] <- "start";
      names[names == "tcnEnd"] <- "end";
      colnames(segs) <- names;
    }

    # Drop bootstrap columns, if any
    names <- colnames(segs);
    keep <- (regexpr("_[0-9]+(|[.][0-9]+)%$", names) == -1);
    segs <- segs[,keep];
  }

  segs;
}, private=TRUE) 



###########################################################################/**
# @RdocMethod writeSegments
#
# @title "Writes the table of segments to file"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{name, tags}{Name and optional tags part of the filename}.
#   \item{path}{The directory where the file will be written.}
#   \item{addHeader}{If @TRUE, header comments are written.}
#   \item{createdBy}{A header comment of whom created the file.}
#   \item{splitters}{If @TRUE, each chromosome is separated by a row
#     of missing values.}
#   \item{overwrite, skip}{If an output file already exists, these
#     arguments specifies what should happen.}
#   \item{...}{Additional arguments pass to \code{getSegments()}.}
# }
#
# \value{
#   Returns the pathname of the the file written.
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
setMethodS3("writeSegments", "PSCBS", function(fit, name=getSampleName(fit), tags=NULL, ext="tsv", path=NULL, addHeader=TRUE, createdBy=NULL, sep="\t", nbrOfDecimals=4L, splitters=FALSE, overwrite=FALSE, skip=FALSE, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'name' and 'tags':
  name <- Arguments$getCharacter(name);
  tags <- Arguments$getCharacters(tags);

  # Argument 'ext':
  ext <- Arguments$getCharacter(ext);

  # Arguments 'path':
  path <- Arguments$getWritablePath(path);

  # Argument 'nbrOfDecimals':
  nbrOfDecimals <- Arguments$getInteger(nbrOfDecimals);


  fullname <- paste(c(name, tags), collapse=",");
  filename <- sprintf("%s.%s", fullname, ext);
  pathname <- Arguments$getWritablePathname(filename, path=path, mustNotExist=(!overwrite && !skip));

  # File already exists?
  if (isFile(pathname)) {
    # Skip?
    if (skip) {
      return(pathname);
    }

    # Overwrite!
    file.remove(pathname);
  }

  # Write to temporary file
  pathnameT <- pushTemporaryFile(pathname);


  sampleName <- getSampleName(fit);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  data <- getSegments(fit, ..., splitters=splitters);

  # Round of floating points
  if (!is.null(nbrOfDecimals)) {
    cols <- tolower(colnames(data));
    isInt <- (regexpr("chromosome|start|end|nbrofloci", cols) != -1);
    cols <- which(isInt);
    for (cc in cols) {
      values <- data[[cc]];
      if (is.double(values)) {
        values <- round(values, digits=0);
        data[[cc]] <- values;
      }
    } # for (key ...)

    cols <- tolower(colnames(data));
    isInt <- (regexpr("chromosome|start|end|nbrofloci", cols) != -1);
    isLog <- (regexpr("call", cols) != -1);
    isDbl <- (!isInt & !isLog);
    cols <- which(isDbl);
    for (kk in cols) {
      values <- data[[kk]];
      if (is.double(values)) {
        values <- round(values, digits=nbrOfDecimals);
        data[[kk]] <- values;
      }
    } # for (key ...) 
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Build header
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (addHeader) {
#    sigmaDelta <- estimateStandardDeviation(fit, method="diff");
   sigmaDelta <- NA;
#    sigmaResiduals <- estimateStandardDeviation(fit, method="res");

    createdOn <- format(Sys.time(), format="%Y-%m-%d %H:%M:%S %Z");
    hdr <- c(
      name=name,
      tags=tags,
      fullname=fullname,
      segmentationMethod=sprintf("segment() of %s", attr(fit, "pkgDetails")),
      nbrOfLoci=nbrOfLoci(fit),
      nbrOfSegments=nbrOfSegments(fit),
      joinSegments=fit$params$joinSegments,
#      signalType=getSignalType(fit),
      sigmaDelta=sprintf("%.4f", sigmaDelta),
#      sigmaResiduals=sprintf("%.4f", sigmaResiduals),
      createdBy=createdBy,
      createdOn=createdOn,
      nbrOfDecimals=nbrOfDecimals,
      nbrOfColumns=ncol(data),
      columnNames=paste(colnames(data), collapse=", "),
      columnClasses=paste(sapply(data, FUN=function(x) class(x)[1]), collapse=", ")
    );
    bfr <- paste("# ", names(hdr), ": ", hdr, sep="");

    cat(file=pathnameT, bfr, sep="\n");
  } # if (addHeader)

  write.table(file=pathnameT, data, append=TRUE, quote=FALSE, sep=sep, 
                                          row.names=FALSE, col.names=TRUE);

  pathname <- popTemporaryFile(pathnameT);

  pathname;  
}) # writeSegments()



############################################################################
# HISTORY:
# 2011-12-03
# o Added writeSegments() for PairedPSCBS.
# o Added argument 'simplify' to getSegments().
# 2011-10-16
# o Added isSegmentSplitter().
# 2011-10-02
# o Now the CBS class extends the AbstractCBS class.
# o Added print() and as.data.frame() to PSCBS.
# o Added getSegments() to PSCBS.
# o DOCUMENTATION: Added Rdoc for several PSCBS methods.
# o Added a PSCBS constructor with documentation.
# 2010-12-01
# o Now also extractByChromosomes() and append() for PSCBS recognizes
#   fields 'tcnLociToExclude' and 'dhLociToExclude'.
# o BUG FIX: extractByChromosome() for PSCBS would call it self instead
#   of extractByChromosomes().
# 2010-11-26
# o Added extractByChromosomes() for PSCBS.
# 2010-09-26
# o getChromosomes() no longer returns NA divers.
# 2010-09-24
# o Added append() and more for PSCBS objects.
############################################################################
