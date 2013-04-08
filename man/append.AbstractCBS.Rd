%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AbstractCBS.RESTRUCT.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 \name{append.AbstractCBS}
\alias{append.AbstractCBS}
\alias{AbstractCBS.append}
\alias{append,AbstractCBS-method}

 \title{Appends one segmentation result to another}

 \description{
   Appends one segmentation result to another,
   where both holds segmentation results \emph{of the same sample}.
 }

 \usage{
\method{append}{AbstractCBS}(...)
}

 \arguments{
  \item{x, other}{The two \code{\link{AbstractCBS}} objects to be combined.}
  \item{addSplit}{If \code{\link[base:logical]{TRUE}}, a "divider" is added between chromosomes.}
  \item{...}{Not used.}
 }

 \value{
   Returns a object of the same class as argument \code{x}.
 }

 \author{Henrik Bengtsson}

 \seealso{
   For more information see \code{\link{AbstractCBS}}.
 }
\keyword{internal}
\keyword{methods}
