%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  PairedPSCBS.EXTS3.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 \name{dropByRegions.PairedPSCBS}
\alias{dropByRegions.PairedPSCBS}
\alias{PairedPSCBS.dropByRegions}
\alias{dropByRegions.PairedPSCBS}
\alias{dropByRegions,PairedPSCBS-method}


 \title{Drops chromosomal regions (a connected set of segments)}

 \description{
   Drops chromosomal regions (a connected set of segments) each of a cetain size (number of segments).
 }

 \usage{\method{dropByRegions}{PairedPSCBS}(this, regions, H=1, ..., verbose=FALSE)}

 \arguments{
  \item{regions}{An \code{\link[base]{integer}} \code{\link[base]{vector}} of length R specifying the indices
    of the left most segment in each of the R regions to be dropped.}
  \item{H}{A non-negative \code{\link[base]{integer}} specifying the size of each region,
    i.e. the number of segments per region.}
  \item{...}{Additional arguments passed to \code{\link[PSCBS:extractByRegions.PairedPSCBS]{*extractByRegions}()}.}
  \item{verbose}{A \code{\link[base]{logical}} or a \code{\link[R.utils]{Verbose}} object.}
 }

 \value{
   Returns a \code{\link{PairedPSCBS}} with (at most) R*H segments dropped.
   If some regions overlap (share segments), then fewer than R*H segments
   are dropped.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   For more information see \code{\link{PairedPSCBS}}.
 }
\keyword{internal}
\keyword{methods}
