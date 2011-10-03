%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  PSCBS.EXTS.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 \name{append.PSCBS}
\alias{append.PSCBS}
\alias{PSCBS.append}
\alias{append.PSCBS}
\alias{append,PSCBS-method}


 \title{Appends one segmentation result to another}

 \description{
   Appends one segmentation result to another.
 }

 \usage{\method{append}{PSCBS}(x, other, addSplit=TRUE, ...)}

 \arguments{
  \item{x, other}{The two \code{\link{PSCBS}} objects to be combined.}
  \item{other}{A \code{\link{PSCBS}} object.}
  \item{addSplit}{If \code{\link[base:logical]{TRUE}}, a "divider" is added between chromosomes.}
  \item{...}{Not used.}
 }

 \value{
   Returns a \code{\link{PSCBS}} object of the same class as argument \code{x}.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   For more information see \code{\link{PSCBS}}.
 }
\keyword{internal}
\keyword{methods}
