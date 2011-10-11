setMethodS3("writeLocusData", "CBS", function(fit, filename=sprintf("%s,byLocus.tsv", getSampleName(fit)), path=NULL, sep="\t", nbrOfDecimals=4L, addHeader=TRUE, createdBy=NULL, overwrite=FALSE, skip=FALSE, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Arguments 'filename' & 'path':
  pathname <- Arguments$getWritablePathname(filename, path=path, mustNotExist=(!overwrite && !skip));

  # Argument 'nbrOfDecimals':
  nbrOfDecimals <- Arguments$getInteger(nbrOfDecimals);


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


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  data <- getLocusData(fit, ...);

  # Round of floating points
  if (!is.null(nbrOfDecimals)) {
    cols <- colnames(data);
    for (key in cols) {
      values <- data[[key]];
      if (is.double(values)) {
        values <- round(values, digits=nbrOfDecimals);
        data[[key]] <- values;
      }
    } # for (key ...)
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Build header
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (addHeader) {
    sigmaDelta <- estimateStandardDeviation(fit, method="diff");
#    sigmaResiduals <- estimateStandardDeviation(fit, method="res");

    createdOn <- format(Sys.time(), format="%Y-%m-%d %H:%M:%S %Z");
    hdr <- c(
      sampleName=getSampleName(fit),
      segmentationMethod=sprintf("segment() of %s", attr(fit, "pkgDetails")),
      nbrOfLoci=nbrOfLoci(fit),
      nbrOfSegments=nbrOfSegments(fit),
      joinSegments=fit$params$joinSegments,
      signalType=getSignalType(fit),
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
}, protected=TRUE) # writeLocusData()



###########################################################################/**
# @set "class=CBS"
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
#   \item{filename, path}{The filename and the path of the file to be written.}
#   \item{addHeader}{If @TRUE, header comments are written.}
#   \item{createdBy}{A header comment of whom created the file.}
#   \item{splitters}{If @TRUE, each chromosome is separated by a row
#     of missing values.}
#   \item{overwrite, skip}{If an output file already exists, these
#     arguments specifies what should happen.}
#   \item{...}{Not used.}
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
setMethodS3("writeSegments", "CBS", function(fit, filename=sprintf("%s.tsv", getSampleName(fit)), path=NULL, addHeader=TRUE, createdBy=NULL, sep="\t", nbrOfDecimals=4L, splitters=FALSE, overwrite=FALSE, skip=FALSE, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Arguments 'filename' & 'path':
  pathname <- Arguments$getWritablePathname(filename, path=path, mustNotExist=(!overwrite && !skip));

  # Argument 'nbrOfDecimals':
  nbrOfDecimals <- Arguments$getInteger(nbrOfDecimals);


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
  data <- getSegments(fit, splitters=splitters);

  # Round of floating points
  if (!is.null(nbrOfDecimals)) {
    cols <- c("start", "end");
    for (key in cols) {
      values <- data[[key]];
      if (is.double(values)) {
        values <- round(values, digits=0);
        data[[key]] <- values;
      }
    } # for (key ...)

    cols <- colnames(data);
    cols <- setdiff(cols, c("chromosome", "start", "end", "nbrOfLoci"));
    for (key in cols) {
      values <- data[[key]];
      if (is.double(values)) {
        values <- round(values, digits=nbrOfDecimals);
        data[[key]] <- values;
      }
    } # for (key ...) 
 }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Build header
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (addHeader) {
    sigmaDelta <- estimateStandardDeviation(fit, method="diff");
#    sigmaResiduals <- estimateStandardDeviation(fit, method="res");

    createdOn <- format(Sys.time(), format="%Y-%m-%d %H:%M:%S %Z");
    hdr <- c(
      sampleName=sampleName,
      segmentationMethod=sprintf("segment() of %s", attr(fit, "pkgDetails")),
      nbrOfLoci=nbrOfLoci(fit),
      nbrOfSegments=nbrOfSegments(fit),
      joinSegments=fit$params$joinSegments,
      signalType=getSignalType(fit),
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
# 2011-09-04
# o Added writeSegments() for CBS.
# o Added writeLocusData() for CBS.
############################################################################
