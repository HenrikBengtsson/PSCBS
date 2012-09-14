%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AbstractCBS.PRUNE.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{pruneByDP.AbstractCBS}
\alias{pruneByDP.AbstractCBS}
\alias{AbstractCBS.pruneByDP}
\alias{pruneByDP.AbstractCBS}
\alias{pruneByDP,AbstractCBS-method}


\title{Prunes the CN profile using dynamical programming}

\description{
 Prunes the CN profile using dynamical programming by specifying the target number of segments or alternative
 how of many change points to drop.
}

\usage{\method{pruneByDP}{AbstractCBS}(fit, nbrOfSegments, ..., verbose=FALSE)}

\arguments{
 \item{nbrOfSegments}{An \code{\link[base]{integer}} specifying the number of segments after
    pruning. If negative, the it specifies the number of change points
    to drop.}
 \item{...}{Optional arguments passed to \code{\link[PSCBS:seqOfSegmentsByDP.AbstractCBS]{*seqOfSegmentsByDP}()}.}
 \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
  Returns a pruned object of the same class.
}

\examples{\dontrun{
 # Drop two segments
 fitP <- pruneByDP(fit, nbrOfSegments=-2);
}}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\references{
  [1] ... \cr
}


\keyword{internal}
\keyword{methods}
