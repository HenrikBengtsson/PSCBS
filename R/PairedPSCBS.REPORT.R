###########################################################################/**
# @set class=PairedPSCBS
# @RdocMethod report
#
# @title "Generates a report of the Paired PSCBS results"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{fit}{A @see "PairedPSCBS" object as returned by 
#      @see "segmentByPairedPSCBS".}
#   \item{sampleName}{A @character string specifying the name of the 
#      tumor-normal pair segmented.}
#   \item{studyName}{A @character string specifying the name of study/project.}
#   \item{...}{Optional arguments passed to the RSP template.}
#   \item{rootPath}{The root directory where to write the report.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   Returns (invisibly) the pathname of the generated PDF.
# }
#
# @examples "../incl/segmentByPairedPSCBS,report.Rex"
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @keyword internal
#*/###########################################################################
setMethodS3("report", "PairedPSCBS", function(fit, sampleName=getSampleName(fit), studyName, ..., rootPath="reports/", verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'sampleName':
  sampleName <- Arguments$getCharacter(sampleName);
  stopifnot(!is.na(sampleName));

  # Argument 'studyName':
  studyName <- Arguments$getCharacter(studyName);

  # Argument 'rootPath':
  rootPath <- Arguments$getWritablePath(rootPath);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Generating Paired PSCBS report");

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

  filenames <- c(
    rsp="PairedPSCBSReport.tex.rsp",
    "PSCBS.bib", "bioinformatics-journals-abbr.bib",
    "natbib.bst"
  );
  verbose && cat(verbose, "Template files:");
  verbose && print(verbose, filenames);

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

  invisible(pdfPathname);
}, protected=TRUE)



############################################################################
# HISTORY:
# 2012-02-27
# o Added Rdoc help.
# o Added report() for PairedPSCBS.
# o Created.
############################################################################
