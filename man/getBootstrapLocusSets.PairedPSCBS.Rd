%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  PairedPSCBS.BOOT.sets.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{getBootstrapLocusSets.PairedPSCBS}
\alias{getBootstrapLocusSets.PairedPSCBS}
\alias{PairedPSCBS.getBootstrapLocusSets}
\alias{getBootstrapLocusSets,PairedPSCBS-method}
\alias{getBootstrapLocusSets}

\title{Generates original and bootstrapped segment-specific index sets}

\description{
 Generates original and bootstrapped segment-specific index sets, which can be used to calculate various bootstrap summaries,
 e.g. segment mean levels.
}

\usage{
\method{getBootstrapLocusSets}{PairedPSCBS}(fit, B=1000L, by=c("betaTN", "betaT"), seed=NULL, verbose=FALSE,
  .validate=FALSE, ...)
}

\arguments{
  \item{B}{A non-negative \code{\link[base]{integer}} specifying the number of bootstrap samples.}
  \item{by}{Should \code{betaTN} or \code{betaT} be used?}
  \item{seed}{An (optional) \code{\link[base]{integer}} specifying the random seed to be
    set before sampling indices.  The random seed is set to its original
    state when exiting.  If \code{\link[base]{NULL}}, it is not set.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
  \item{.validate}{If \code{\link[base:logical]{TRUE}}, additional sanity checks are performed
    to validate the correctness.  This is only needed for troubleshooting
    if it is suspected there is a bug.}
  \item{...}{Not used.}
}

\value{
  Returns a \code{\link[base]{list}}.
}

\author{Henrik Bengtsson}

\seealso{
  This is used internally by various bootstrap methods.
}


\keyword{internal}
\keyword{methods}
