%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  PairedPSCBS.REPORT.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{report.PairedPSCBS}
\alias{report.PairedPSCBS}
\alias{PairedPSCBS.report}
\alias{report.PairedPSCBS}
\alias{report,PairedPSCBS-method}


\title{Generates a report of the Paired PSCBS results}

\description{
 Generates a report of the Paired PSCBS results.
}

\usage{\method{report}{PairedPSCBS}(fit, sampleName=getSampleName(fit), studyName, ..., rootPath="reports/", verbose=FALSE)}

\arguments{
  \item{fit}{A \code{\link{PairedPSCBS}} object as returned by
     \code{\link{segmentByPairedPSCBS}}().}
  \item{sampleName}{A \code{\link[base]{character}} string specifying the name of the
     tumor-normal pair segmented.}
  \item{studyName}{A \code{\link[base]{character}} string specifying the name of study/project.}
  \item{...}{Optional arguments passed to the RSP template.}
  \item{rootPath}{The root directory where to write the report.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
  Returns (invisibly) the pathname of the generated PDF.
}

\examples{

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Load SNP microarray data
# (note to package developers: this example data set may
#  be replaced in a future release of the package)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
pathname <- system.file("data-ex/PairedPSCBS,exData,chr01.Rbin", package="PSCBS")
data <- R.utils::loadObject(pathname)
str(data)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Paired PSCBS segmentation
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Drop single-locus outliers
dataS <- dropSegmentationOutliers(data)

# Speed up example by segmenting fewer loci
dataS <- dataS[seq(from=1, to=nrow(data), by=5),]

str(dataS)

R.oo::attachLocally(dataS)

# Paired PSCBS segmentation
fit <- segmentByPairedPSCBS(CT, betaT=betaT, betaN=betaN,
                            chromosome=chromosome, x=x,
                            seed=0xBEEF, verbose=-10)

report(fit, sampleName="PairedPSCBS", studyName="PSCBS-Ex", verbose=-10)
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{PairedPSCBS}}.
}


\keyword{internal}
\keyword{methods}
