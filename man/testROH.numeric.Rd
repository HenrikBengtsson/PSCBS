%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  testROH.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{testROH.numeric}
\alias{testROH.numeric}

\title{Tests if a segment is in Run-of-Homozygosity (ROH)}

\description{
 Tests if a segment is in Run-of-Homozygosity (ROH).
}

\usage{
\method{testROH}{numeric}(betaN, muN, csN=NULL, minNbrOfSnps=1, delta=1/12, ..., verbose=FALSE)
}

\arguments{
  \item{betaN}{A \code{\link[base]{numeric}} \code{\link[base]{vector}} of J matched normal BAFs in [0,1]
       (due to noise, values may be slightly outside as well) or \code{\link[base]{NA}}
       for non-polymorphic loci.}
  \item{muN}{An optional \code{\link[base]{numeric}} \code{\link[base]{vector}} of J genotype calls in
       \{0,1/2,1\} for AA, AB, and BB, respectively,
       and \code{\link[base]{NA}} for non-polymorphic loci.}
  \item{csN}{A \code{\link[base]{numeric}} \code{\link[base]{vector}} of J genotype confidence scores.
       If \code{\link[base]{NULL}}, ad hoc scores calculated from \code{betaN} are used.}
  \item{minNbrOfSnps}{Minimum number of SNPs required to test segment.
       If not tested, \code{\link[base]{NA}} is returned.}
  \item{delta}{A \code{\link[base]{double}} scalar specifying the maximum (weighted)
       proportion of heterozygous SNPs allowed in an ROH region.}
  \item{...}{Not used.}
 \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
 Returns a \code{\link[base]{logical}}.
}

\author{Pierre Neuvial, Henrik Bengtsson}


\keyword{methods}
\keyword{internal}
