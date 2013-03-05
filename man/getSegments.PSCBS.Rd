%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  PSCBS.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{getSegments.PSCBS}
\alias{getSegments.PSCBS}
\alias{PSCBS.getSegments}
\alias{getSegments,PSCBS-method}

 \title{Gets the segments}

 \description{
   Gets the segments.
 }

 \usage{
\method{getSegments}{PSCBS}(fit, simplify=FALSE, splitters=TRUE, addGaps=FALSE, ...)
}

 \arguments{
  \item{simplify}{If \code{\link[base:logical]{TRUE}}, redundant and intermediate information is dropped.}#  \item{splitters}{If \code{\link[base:logical]{TRUE}}, "splitters" between chromosomes are
     preserved, otherwise dropped.}
  \item{...}{Not used.}
 }

 \value{
   Returns a SxK \code{\link[base]{data.frame}}, where S in the number of segments,
   and K is the number of segment-specific fields.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   For more information see \code{\link{PSCBS}}.
 }
\keyword{internal}
\keyword{methods}
