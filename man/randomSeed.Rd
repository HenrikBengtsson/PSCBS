%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  randomSeed.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{randomSeed}
\alias{randomSeed}


\title{Sets and resets the .Random.seed in the global environment}

\description{
 Sets and resets the .Random.seed in the global environment.
}

\usage{
 randomSeed(action=c("set", "reset", "get"), seed=NULL)
}

\arguments{
  \item{action}{A \code{\link[base]{character}} string specifying the action.}
  \item{seed}{Random seed to be set; only for \code{action="set"}.}
}

\value{
  Returns invisibly the previous/current \code{.Random.seed}.
}

\author{Henrik Bengtsson}


\keyword{IO}
