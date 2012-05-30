###########################################################################/**
# @set class=AbstractCBS
# @RdocMethod report
#
# @title "Generates a report of the segmentation results"
#
# \description{
#  @get "title".
#  Currently reports can be generated for segmentation results of class
#  @see "CBS" and @see "PairedPSCBS".
# }
#
# @synopsis
#
# \arguments{
#   \item{fit}{An @see "AbstractCBS" object.}
#   \item{sampleName}{A @character string specifying the name of the 
#      sample segmented.}
#   \item{studyName}{A @character string specifying the name of study/project.}
#   \item{...}{Optional arguments passed to the RSP template.}
#   \item{rootPath}{The root directory where to write the report.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   Returns the pathname of the generated PDF.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @keyword internal
#*/###########################################################################
setMethodS3("report", "AbstractCBS", function(fit, sampleName=getSampleName(fit), studyName, ..., rootPath="reports/", .filenames=c(rsp="*", "PSCBS.bib", "bioinformatics-journals-abbr.bib", "natbib.bst"), verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'sampleName':
  sampleName <- Arguments$getCharacter(sampleName);
  if (is.na(sampleName)) {
    throw("Cannot generate report. Argument 'sampleName' is non-valid or missing.");
  }

  # Argument 'studyName':
  if (is.missing(studyName)) {
    throw("Cannot generate report. Argument 'studyName' is missing.");
  }
  studyName <- Arguments$getCharacter(studyName);
  if (is.na(studyName)) {
    throw("Cannot generate report. Argument 'studyName' is non-valid.");
  }

  # Argument 'rootPath':
  rootPath <- Arguments$getWritablePath(rootPath);

  # Argument '.filenames':
  if (!is.null(.filenames)) {
    .filenames <- Arguments$getCharacters(.filenames, useNames=TRUE);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Generating CBS report");

  verbose && cat(verbose, "Sample name: ", sampleName);
  verbose && cat(verbose, "Number of chromosomes: ", nbrOfChromosomes(fit));
  verbose && cat(verbose, "Number of segments: ", nbrOfSegments(fit));

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Report template arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Default arguments
  rspArgs <- list(
    fit = fit,
    sampleName = sampleName,
    studyName = studyName,
    dataSet = NULL,
    Clim = c(0,4),
    Blim = c(0,1),
    figForce = FALSE
  );

  # Override with user arguments
  userArgs <- list(...);
  for (key in names(userArgs)) {
    rspArgs[[key]] <- userArgs[[key]];
  }

  if (is.null(rspArgs$reportPath)) {
    rspArgs$reportPath <- file.path(rootPath, rspArgs$studyName);
  }
  rspArgs$reportPath <- Arguments$getWritablePath(rspArgs$reportPath);
  verbose && cat(verbose, "Report root path: ", rspArgs$reportPath);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Copy report files
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Copy report files");

  # Directory where all report templates lives
  srcPath <- system.file("templates", "rsp", package="PSCBS");
  srcPath <- Arguments$getReadablePath(srcPath);
  verbose && cat(verbose, "Source path: ", srcPath);

  filenames <- .filenames;
  idx <- which(filenames == "*");
  if (length(idx) > 0) {
    className <- class(fit)[1];
    rsp <- sprintf("%sReport.tex.rsp", className);
    filenames[idx] <- rsp;
  }
  verbose && cat(verbose, "Template files:");
  verbose && print(verbose, filenames);

  # Sanity check
  stopifnot(is.element("rsp", names(filenames)));

  destFilenames <- filenames;
  destFilenames["rsp"] <- sprintf("%s,%s", sampleName, filenames["rsp"]);

  for (kk in seq(along=filenames)) {
    filename <- filenames[kk];
    destFilename <- destFilenames[kk];
    verbose && enter(verbose, sprintf("File #%d ('%s') of %d", kk, filename, length(filenames)));
    srcPathname <- Arguments$getReadablePathname(filename, path=srcPath);
    pathname <- filePath(rspArgs$reportPath, destFilename);
    verbose && cat(verbose, "Source: ", srcPathname);
    verbose && cat(verbose, "Destination: ", pathname);
    if (!isFile(pathname)) {
      copyFile(srcPathname, pathname, ...);
      # Sanity check
      stopifnot(isFile(pathname));
    }
    verbose && exit(verbose);
  }

  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Build reports
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Processing RSP template");
  opwd <- setwd(rspArgs$reportPath);
  on.exit({
    if (!is.null(opwd)) setwd(opwd);
  }, add=TRUE);
  rspArgs$reportPath <- ".";
  rspArgs$figPath <- "figures/";

  filename <- destFilenames["rsp"];
  pdfPathname <- R.rsp::rsp(filename, ..., verbose=verbose);
  pdfPathname <- getAbsolutePath(pdfPathname);
  verbose && exit(verbose);

  setwd(opwd);
  opwd <- NULL;
  pdfPathname <- getRelativePath(pdfPathname);
  verbose && cat(verbose, "Final report: ", pdfPathname);

  verbose && exit(verbose);

  pdfPathname;
}, protected=TRUE)



############################################################################
# HISTORY:
# 2012-05-30
# o Now report() gives more a informative error message if arguments
#   'sampleName' or 'studyName' are non-valid or missing.
# 2012-04-20
# o Added argument '.filenames'.
# o Created from former PairedPSCBS.REPORT.R, which history as below.
# 2012-02-27
# o Added Rdoc help.
# o Added report() for PairedPSCBS.
# o Created.
############################################################################
