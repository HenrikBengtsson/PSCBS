%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  CBS.RESTRUCT.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 \name{append.CBS}
\alias{append.CBS}
\alias{CBS.append}
\alias{append,CBS-method}

 \title{Appends one segmentation result to another}

 \description{
   Appends one segmentation result to another.
 }

 \usage{
\method{append}{CBS}(x, other, addSplit=TRUE, ...)
}

 \arguments{
  \item{x, other}{The two \code{\link{CBS}} objects to be combined.}
  \item{other}{A \code{\link{PSCBS}} object.}
  \item{addSplit}{If \code{\link[base:logical]{TRUE}}, a "divider" is added between chromosomes.}
  \item{...}{Not used.}
 }

 \value{
   Returns a \code{\link{CBS}} object of the same class as argument \code{x}.
 }

 \author{Henrik Bengtsson}

 \seealso{
   For more information see \code{\link{CBS}}.
 }
\keyword{internal}
\keyword{methods}
