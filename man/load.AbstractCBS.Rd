%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AbstractCBS.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{AbstractCBS$load}
\alias{AbstractCBS$load}
\alias{load.AbstractCBS}
\alias{AbstractCBS.load}
\alias{load,AbstractCBS-method}
\alias{load}

\title{Loads an AbstractCBS object from file}

\description{
 Loads an AbstractCBS object from file and assert that it is of the requested class.
}

\usage{
## Static method (use this):
## AbstractCBS$load(...)

## Don't use the below:
\method{load}{AbstractCBS}(static, ...)
}

\arguments{
  \item{...}{Additional arguments passed to \code{\link[R.utils]{loadObject}}.}
}

\value{
  Returns the loaded AbstractCBS object.
}

\author{Henrik Bengtsson}

\seealso{
  Internally \code{\link[R.utils]{loadObject}} is used.
  To save an object, see \code{\link[PSCBS:save.AbstractCBS]{*save}()}.
  For more information see \code{\link{AbstractCBS}}..
}
\keyword{internal}
\keyword{methods}
