%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AbstractCBS.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{resetSegments.AbstractCBS}
\alias{resetSegments.AbstractCBS}
\alias{AbstractCBS.resetSegments}
\alias{resetSegments,AbstractCBS-method}

 \title{Reset the segments}

 \description{
   Reset the segments.  More precisely, it removes columns in the segmentation
   result table that have been added by methods after the actual
   segmentation method, e.g. bootstrap estimated mean level quantiles
   and various calls.
   It leave the basic segmentation results untouched,
   i.e. the partitioning and the segment means.
 }

 \usage{
\method{resetSegments}{AbstractCBS}(fit, ...)
}

 \arguments{
  \item{...}{Not used.}
 }

 \value{
   Returns an object if the same class as the input result.
 }

 \author{Henrik Bengtsson}

 \seealso{
   For more information see \code{\link{AbstractCBS}}.
 }
\keyword{internal}
\keyword{methods}
