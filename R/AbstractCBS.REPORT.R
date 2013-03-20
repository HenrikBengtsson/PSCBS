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
#   \item{rspTags}{Optional @character @vector of tags for further specifying
#      which RSP report to generate.}
#   \item{rootPath}{The root directory where to write the report.}
#   \item{force}{If @TRUE, RSP template files are copied to the reports/
#      directory, regardless if they already exists there or not.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   Returns the pathname of the generated PDF.
# }
#
# @author "HB"
#
# \seealso{
#   @seeclass
# }
#
# @keyword internal
#*/###########################################################################
setMethodS3("report", "AbstractCBS", function(fit, sampleName=getSampleName(fit), studyName, ..., rspTags=NULL, rootPath="reports/", .filenames=c(rsp="*", "PSCBS.bib", "bioinformatics-journals-abbr.bib", "natbib.bst"), force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'sampleName':
  sampleName <- Arguments$getCharacter(sampleName);
  if (is.na(sampleName)) {
    throw("Cannot generate report. Argument 'sampleName' is non-valid or missing.");
  }

  # Argument 'studyName':
  if (missing(studyName)) {
    throw("Cannot generate report. Argument 'studyName' is missing.");
  }
  studyName <- Arguments$getCharacter(studyName);
  if (is.na(studyName)) {
    throw("Cannot generate report. Argument 'studyName' is non-valid.");
  }

  # Argument 'rspTags':
  if (!is.null(rspTags)) {
    rspTags <- Arguments$getCharacters(rspTags);
    rspTags <- unlist(strsplit(rspTags, split=",", fixed=TRUE));
    rspTags <- rspTags[nchar(rspTags) > 0L];
  }

  # Argument 'rootPath':
  rootPath <- Arguments$getWritablePath(rootPath);

  # Argument '.filenames':
  if (!is.null(.filenames)) {
    .filenames <- Arguments$getCharacters(.filenames, useNames=TRUE);
  }

  # Argument 'force':
  force <- Arguments$getLogical(force);

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
  srcPath <- "templates";

  # If missing, default to one that comes with PSCBS/templates/
  if (!isDirectory(srcPath)) {
    srcPath <- system.file("templates", package="PSCBS");
  }
  srcPath <- file.path(srcPath, "rsp");
  srcPath <- Arguments$getReadablePath(srcPath);
  verbose && cat(verbose, "Source path: ", srcPath);

  filenames <- .filenames;

  # Construct the filename of the main RSP file to compile, iff missing.
  idx <- which(filenames == "*");
  if (length(idx) > 0) {
    className <- class(fit)[1];
    fullname <- paste(c(className, rspTags), collapse=",");
    rsp <- sprintf("%s,report.tex.rsp", fullname);
    filenames[idx] <- rsp;
  }

  # Include files that are listed in the optional '.install_extras' file.
  extras <- file.path(srcPath, ".install_extras");
  if (isFile(extras)) {
    extras <- readLines(extras, warn=FALSE);
    extras <- extras[sapply(file.path(srcPath, extras), FUN=isFile)];
    filenames <- c(filenames, extras);
  }

  verbose && cat(verbose, "Template files:");
  verbose && print(verbose, filenames);

  # Sanity check
  stopifnot(is.element("rsp", names(filenames)));

  # Make sure 'rsp' is first
  idx <- which(names(filenames) == "rsp");
  filename <- c(filenames[idx], filenames[-idx]);

  # Drop duplicated filenames
  filenames <- filenames[!duplicated(filenames)];

  rspFilename <- filenames["rsp"];
  rspPathname <- file.path(srcPath, rspFilename);
  verbose && cat(verbose, "RSP template file: ", rspPathname);
  rspPathname <- Arguments$getReadablePathname(rspPathname);

  destFilenames <- filenames;
  destFilenames["rsp"] <- sprintf("%s,%s", sampleName, rspFilename);

  for (kk in seq(along=filenames)) {
    filename <- filenames[kk];
    destFilename <- destFilenames[kk];
    verbose && enter(verbose, sprintf("File #%d ('%s') of %d", kk, filename, length(filenames)));
    srcPathname <- Arguments$getReadablePathname(filename, path=srcPath);
    pathname <- filePath(rspArgs$reportPath, destFilename);
    verbose && cat(verbose, "Source: ", srcPathname);
    verbose && cat(verbose, "Destination: ", pathname);
    if (force || !isFile(pathname)) {
      copyFile(srcPathname, pathname, overwrite=force, ...);
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
# 2013-03-09
# o Now report() also inclused files listed in the optional file
#   '.install_extras' of the source RSP template directory.
#   The same filename is used by 'R CMD build/check' for including
#   extra vignette source files.
# 2012-09-18
# o Added argument 'force' to report() for AbstractCBS.  This will
#   copy the RSP template files again, although they are already in
#   reports/ output directory.
# o Now report(fit, ..., rspTags) for AbstractCBS looks for the RSP 
#   template named <className>(,<rspTags>),report.tex.rsp, where
#   className is class(fit)[1] and  argument 'rspTags' is an optional
#   comma-separated character string/vector.
# o Now report() for AbstractCBS looks for the RSP template in templates/,
#   and as a backup in templates,PSCBS/.  If the latter does not exist,
#   it is automatically created as a soft link to templates/ of the
#   PSCBS package.  This allows anyone to create their own customized
#   copy (in templates/) of the default PSCBS RSP report.
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
