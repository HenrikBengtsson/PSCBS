%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AbstractCBS.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{getSegments.AbstractCBS}
\alias{getSegments.AbstractCBS}
\alias{AbstractCBS.getSegments}
\alias{getSegments.AbstractCBS}
\alias{getSegments,AbstractCBS-method}


 \title{Gets the segments}

 \description{
   Gets the segments.
 }

 \usage{\method{getSegments}{AbstractCBS}(...)}

 \arguments{
  \item{splitters}{If \code{\link[base:logical]{TRUE}}, "splitters" between chromosomes are
     preserved, otherwise dropped.}
  \item{...}{Not used.}
 }

 \value{
   Returns a SxK \code{\link[base]{data.frame}}, where S in the number of segments,
   and K is the number of segment-specific fields.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   For more information see \code{\link{AbstractCBS}}.
 }
\keyword{internal}
\keyword{methods}
