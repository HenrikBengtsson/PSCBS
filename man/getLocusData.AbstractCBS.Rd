%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AbstractCBS.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{getLocusData.AbstractCBS}
\alias{getLocusData.AbstractCBS}
\alias{AbstractCBS.getLocusData}
\alias{getLocusData,AbstractCBS-method}

 \title{Gets the locus-level data}

 \description{
   Gets the locus-level data.
 }

 \usage{
\method{getLocusData}{AbstractCBS}(...)
}

 \arguments{
  \item{splitters}{If \code{\link[base:logical]{TRUE}}, "splitters" between chromosomes are
     preserved, otherwise dropped.}
  \item{...}{Not used.}
 }

 \value{
   Returns a JxL \code{\link[base]{data.frame}}, where J in the number of loci,
   and L is the number of locus-specific fields.
 }

 \author{Henrik Bengtsson}

 \seealso{
   For more information see \code{\link{AbstractCBS}}.
 }
\keyword{internal}
\keyword{methods}
