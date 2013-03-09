%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AbstractCBS.RESTRUCT.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{mergeThreeSegments.AbstractCBS}
\alias{mergeThreeSegments.AbstractCBS}
\alias{AbstractCBS.mergeThreeSegments}
\alias{mergeThreeSegments,AbstractCBS-method}

 \title{Merge a segment and its two flanking segments}

 \description{
   Merge a segment and its two flanking segments into one segment, and recalculating the segment statistics.
 }

 \usage{
\method{mergeThreeSegments}{AbstractCBS}(fit, middle, ...)
}

 \arguments{
  \item{middle}{An \code{\link[base]{integer}} specifying the three segments
    (middle-1, middle, middle+1) to be merged.}
  \item{...}{Additional arguments passed to \code{\link[PSCBS:mergeTwoSegments.AbstractCBS]{*mergeTwoSegments}()}.}
 }

 \value{
   Returns an \code{\link{AbstractCBS}} of the same class with two less segment.
 }

 \author{Henrik Bengtsson}

 \seealso{
   Internally \code{\link[PSCBS:mergeTwoSegments.AbstractCBS]{*mergeTwoSegments}()} is used.
   For more information see \code{\link{AbstractCBS}}.
 }
\keyword{internal}
\keyword{methods}
