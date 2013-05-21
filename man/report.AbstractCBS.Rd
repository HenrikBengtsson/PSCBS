%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AbstractCBS.REPORT.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{report.AbstractCBS}
\alias{report.AbstractCBS}
\alias{AbstractCBS.report}
\alias{report,AbstractCBS-method}

\title{Generates a report of the segmentation results}

\description{
 Generates a report of the segmentation results.
 Currently reports can be generated for segmentation results of class
 \code{\link{CBS}} and \code{\link{PairedPSCBS}}.
}

\usage{
\method{report}{AbstractCBS}(fit, sampleName=getSampleName(fit), studyName, ..., rspTags=NULL,
  rootPath="reports/", .filename="*", envir=new.env(), verbose=FALSE)
}

\arguments{
  \item{fit}{An \code{\link{AbstractCBS}} object.}
  \item{sampleName}{A \code{\link[base]{character}} string specifying the name of the
     sample segmented.}
  \item{studyName}{A \code{\link[base]{character}} string specifying the name of study/project.}
  \item{...}{Optional arguments passed to the RSP template.}
  \item{rspTags}{Optional \code{\link[base]{character}} \code{\link[base]{vector}} of tags for further specifying
     which RSP report to generate.}
  \item{rootPath}{The root directory where to write the report.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
  Returns the pathname of the generated PDF.
}

\author{Henrik Bengtsson}

\seealso{
  For more information see \code{\link{AbstractCBS}}.
}


\keyword{internal}
\keyword{methods}
