%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AbstractCBS.RESTRUCT.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{dropRegions.AbstractCBS}
\alias{dropRegions.AbstractCBS}
\alias{AbstractCBS.dropRegions}
\alias{dropRegions,AbstractCBS-method}
 \alias{AbstractCBS.dropRegion}
\alias{dropRegion.AbstractCBS}
\alias{dropRegion,AbstractCBS-method}


 \title{Drops chromosomal regions (a connected set of segments)}

 \description{
   Drops chromosomal regions (a connected set of segments) each of a certain size (number of segments).
   \emph{None of the statistics are recalculated}.
 }

 \usage{
\method{dropRegions}{AbstractCBS}(this, regions, H=1, ..., asMissing=FALSE, verbose=FALSE)
}

 \arguments{
  \item{regions}{An \code{\link[base]{integer}} \code{\link[base]{vector}} of length R specifying the indices
    of the left most segment in each of the R regions to be dropped.}
  \item{H}{A non-negative \code{\link[base]{integer}} specifying the size of each region,
    i.e. the number of segments per region.}
  \item{...}{Additional arguments passed to \code{\link[PSCBS:extractRegions.AbstractCBS]{*extractRegions}()}.}
  \item{asMissing}{If \code{\link[base:logical]{TRUE}}, dropped segments are replaced by missing values,
    otherwise they are truly dropped.}
  \item{verbose}{A \code{\link[base]{logical}} or a \code{\link[R.utils]{Verbose}} object.}
 }

 \value{
   Returns an \code{\link{AbstractCBS}} object of the same class with (at most)
   R*H segments dropped.
   If some regions overlap (share segments), then fewer than R*H segments
   are dropped.
 }

 \author{Henrik Bengtsson}

 \seealso{
   Internally \code{\link[PSCBS:extractRegions.AbstractCBS]{*extractRegions}()} is used.
   See also \code{\link[PSCBS:dropChangePoint.AbstractCBS]{*dropChangePoint}()} and \code{\link[PSCBS:mergeTwoSegments.AbstractCBS]{*mergeTwoSegments}()}.
   For more information see \code{\link{AbstractCBS}}.
 }
\keyword{internal}
\keyword{methods}
