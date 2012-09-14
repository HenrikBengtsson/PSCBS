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

 \alias{gapsToSegments}

 \title{Gets the genomic segments that are complementary to the gaps}

 \description{
  Gets the genomic segments that are complementary to the gaps, with default chromosome boundaries being \code{-Inf}
  and \code{+Inf}.
 }

 \usage{\method{gapsToSegments}{data.frame}(gaps, resolution=1, minLength=0, dropGaps=FALSE, ...)}

 \arguments{
   \item{gaps}{A \code{\link[base]{data.frame}} with columns \code{chromosome}, \code{start},
     and \code{stop}. Any overlapping gaps will throw an error.}
   \item{resolution}{A non-negative \code{\link[base]{numeric}} specifying the minimum
     length unit, which by default equals one nucleotide/base pair.}
   \item{minLength}{Minimum length of segments to be kept.}
   \item{dropGaps}{If \code{\link[base:logical]{TRUE}}, the gaps themselves are not part of the output.}
   \item{...}{Not used.}
 }

 \value{
   Returns \code{\link[base]{data.frame}} with columns \code{chromosome} (if that argument
   is given), \code{start}, \code{stop} and \code{length}.
   The segments are ordered along the genome.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   \code{\link{findLargeGaps}}().
 }



\keyword{methods}
\keyword{IO}
\keyword{internal}
