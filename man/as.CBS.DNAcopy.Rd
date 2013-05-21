%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  CBS.EXTS.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{as.CBS.DNAcopy}
\alias{as.CBS.DNAcopy}
\alias{DNAcopy.as.CBS}
\alias{as.CBS,DNAcopy-method}

\title{Coerces a DNAcopy object to a CBS object}

\description{
 Coerces a DNAcopy object to a CBS object.
}

\usage{
\method{as.CBS}{DNAcopy}(fit, sample=1L, ...)
}

\arguments{
  \item{fit}{A \code{\link{DNAcopy}} object (of the \pkg{DNAcopy} package.)}
  \item{sample}{An index specifying which sample to extract,
    if more than one exists.}
  \item{...}{Not used.}
}

\value{
  Returns a \code{\link{CBS}} object.
}

\author{Henrik Bengtsson}

\seealso{
  \code{\link[PSCBS:as.DNAcopy.CBS]{as.DNAcopy()}}.
  For more information see \code{\link{DNAcopy}}.
}


\keyword{internal}
\keyword{methods}
