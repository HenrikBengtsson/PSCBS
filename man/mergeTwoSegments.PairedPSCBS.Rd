%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  PairedPSCBS.RESTRUCT.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 \name{mergeTwoSegments.PairedPSCBS}
\alias{mergeTwoSegments.PairedPSCBS}
\alias{PairedPSCBS.mergeTwoSegments}
\alias{mergeTwoSegments,PairedPSCBS-method}

 \title{Merge two neighboring segments}

 \description{
   Merge two neighboring segments by recalculating segment statistics.
 }

 \usage{
\method{mergeTwoSegments}{PairedPSCBS}(this, left, update=TRUE, verbose=FALSE, ...)
}

 \arguments{
  \item{left}{An \code{\link[base]{integer}} specifying the segments (left, left+1)
    to be merged.}
  \item{update}{If \code{\link[base:logical]{TRUE}}, segment statistics are updated.}
  \item{verbose}{A \code{\link[base]{logical}} or a \code{\link[R.utils]{Verbose}} object.}
  \item{...}{Not used.}
 }

 \value{
   Returns a \code{\link{PairedPSCBS}} with one less segment.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   To drop regions (a connected set of segments) see \code{dropRegions()}.
   For more information see \code{\link{PairedPSCBS}}.
 }
\keyword{internal}
\keyword{methods}
