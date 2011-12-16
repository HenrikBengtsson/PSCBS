%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  gapsToSegments.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 \name{gapsToSegments.data.frame}
\alias{gapsToSegments.data.frame}
\alias{gapsToSegments.data.frame}


 \title{Gets the genomic segments that are complementary to the gaps}

 \description{
  Gets the genomic segments that are complementary to the gaps, with default chromosome boundaries being \code{-Inf}
  and \code{+Inf}.
 }

 \usage{\method{gapsToSegments}{data.frame}(gaps, resolution=1, ...)}

 \arguments{
   \item{gaps}{A \code{\link[base]{data.frame}} with columns \code{chromosome}, \code{start},
     and \code{stop}. Any overlapping gaps will throw an error.}
   \item{resolution}{A non-negative \code{\link[base]{numeric}} specifying the minimum
     length unit, which by default equals one nucleotide/base pair.}
   \item{...}{Not used.}
 }

 \value{
   Returns \code{\link[base]{data.frame}} with columns \code{chromosome} (if that argument
   is given), \code{start}, \code{stop} and \code{length}.
   The segments are ordered along the genome.
   Zero length segments are automatically dropped.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   \code{\link{findLargeGaps}}().
 }



\keyword{methods}
\keyword{IO}
\keyword{internal}
