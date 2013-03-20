%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  findCopyNumberNeutralState.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{findNeutralCopyNumberState}
\alias{findNeutralCopyNumberState.default}
\alias{findNeutralCopyNumberState}

\title{Call segments to be copy neutral based on allelic imbalance calls and total copy number estimates}

\description{
 Call segments to be copy neutral based on allelic imbalance calls and total copy number estimates.
}

\usage{
\method{findNeutralCopyNumberState}{default}(C, isAI, weights=NULL, ..., minDensity=1e-10, flavor=c("firstPeak", "maxPeak"), verbose=FALSE)
}

\arguments{
  \item{C}{A \code{\link[base]{numeric}} \code{\link[base]{vector}} of region-level total copy number estimates.}
  \item{isAI}{A \code{\link[base]{logical}} \code{\link[base]{vector}} of "allelic imbalance" calls.}
  \item{weights}{An optional \code{\link[base]{numeric}} \code{\link[base]{vector}} of non-negative weights.}
  \item{...}{Further argumants to be passed to the density estimation
    function.}
  \item{minDensity}{A \code{\link[base]{numeric}} value, below which density peaks are
    discarded.}
  \item{flavor}{A \code{\link[base]{character}} string specifying how to identify the
    mode of the AB segments.}
  \item{verbose}{If \code{\link[base:logical]{TRUE}}, extra information is output.}
}

\value{
  A \code{\link[base]{logical}} \code{\link[base]{vector}} of "neutral copy number state" calls.
}

\author{Pierre Neuvial, Henrik Bengtsson}


\keyword{internal}
