%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  PairedPSCBS.callROH.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{callROH.PairedPSCBS}
\alias{callROH.PairedPSCBS}
\alias{PairedPSCBS.callROH}
\alias{callROH,PairedPSCBS-method}
\alias{callROH.NonPairedPSCBS}

\title{Calls segments that are in ROH}

\description{
 Calls segments that are in ROH, i.e. that have no (true) heterozygous genotypes.
 Run of homozygosity (ROH) is a property of the normal (germline) sample.
}

\usage{
\method{callROH}{PairedPSCBS}(fit, ..., updateMeans=TRUE, force=FALSE, verbose=FALSE)
}

\arguments{
  \item{...}{Additional arguments passed to \code{\link{testROH}}().}
  \item{updateMeans}{If \code{\link[base:logical]{TRUE}}, DH and (C1,C2) mean levels are set
   to \code{\link[base]{NA}} for segments called ROH, otherwise not.}
  \item{force}{If \code{\link[base:logical]{FALSE}}, and ROH calls already exits,
   then nothing is done, otherwise the calls are done.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
  Returns a \code{\link{PairedPSCBS}} object with ROH calls.
}

\author{Pierre Neuvial, Henrik Bengtsson}

\seealso{
  Internally, \code{\link{testROH}}() is used.
  To call allelic balance (AB) see \code{\link[PSCBS:callAB.PairedPSCBS]{*callAB}()}.
  To call loss of heterozygosity (LOH) see \code{\link[PSCBS:callLOH.PairedPSCBS]{*callLOH}()}.
}
\keyword{internal}
\keyword{methods}
