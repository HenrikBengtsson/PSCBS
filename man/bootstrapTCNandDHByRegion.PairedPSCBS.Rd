%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  PairedPSCBS.BOOT2.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{bootstrapTCNandDHByRegion.PairedPSCBS}
\alias{bootstrapTCNandDHByRegion.PairedPSCBS}
\alias{PairedPSCBS.bootstrapTCNandDHByRegion}
\alias{bootstrapTCNandDHByRegion,PairedPSCBS-method}

\title{Estimate confidence intervals of TCN and DH segment levels}

\description{
 Estimate confidence intervals of TCN and DH segment levels using bootstrap.
}

\usage{
\method{bootstrapTCNandDHByRegion}{PairedPSCBS}(fit, B=1000, statsFcn=function(x) quantile(x, probs = c(0.025, 0.05, 0.95, 0.975), na.rm = TRUE), by=c("betaTN", "betaT"), force=FALSE, seed=NULL, verbose=FALSE, .debug=FALSE, ...)
}

\arguments{
  \item{B}{A postive \code{\link[base]{integer}} specifying the number of bootstrap samples.}
  \item{statsFcn}{A \code{\link[base]{function}} that estimates confidence intervals given
     locus-level data.}
  \item{by}{A \code{\link[base]{character}} specifying whether DH should be calculated from
     normalized ('betaTN') or non-normalized ('betaT') tumor BAFs.}
  \item{force}{If \code{\link[base:logical]{TRUE}}, already existing estimates are ignored,
     otherwise not.}
  \item{seed}{(optional) A random seed.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
  \item{.debug}{(internal) If \code{\link[base:logical]{TRUE}}, additional sanity checks are
     performed internally..}
  \item{...}{Not used.}
}

\value{
  Returns a \code{\link{PairedPSCBS}} object with copy-neutral calls.
}

\author{Henrik Bengtsson}

\seealso{
  Internally, one of the following methods are used:
  \code{\link[PSCBS:callCopyNeutralByTCNofAB.PairedPSCBS]{*callCopyNeutralByTCNofAB}()}.
}
\keyword{internal}
\keyword{methods}
