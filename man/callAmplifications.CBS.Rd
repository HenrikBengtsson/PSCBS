%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  CBS.CALL.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{callAmplifications.CBS}
\alias{callAmplifications.CBS}
\alias{CBS.callAmplifications}
\alias{callAmplifications,CBS-method}

\title{Calls (focal) amplifications}

\description{
 Calls (focal) amplifications.
}

\usage{
\method{callAmplifications}{CBS}(fit, adjust=1, maxLength=2e+07, method=c("ucsf-exp"), ...,
  verbose=FALSE)
}

\arguments{
 \item{adjust}{A positive scale factor adjusting the sensitivity of the
   caller, where a value less (greater) than 1.0 makes the caller
   less (more) sensitive.}
 \item{maxLength}{A \code{\link[base]{double}} scalar specifying the maximum length of a segment
   in order for it to be considered a focal amplification.}
 \item{method}{A \code{\link[base]{character}} string specifying the calling algorithm to use.}
 \item{...}{Additional/optional arguments used to override the default
   parameters used by the caller.}
 \item{verbose}{\code{\link[R.utils]{Verbose}}.}
}

\value{
 Returns a \code{\link[PSCBS]{CBS}} object where \code{\link[base]{logical}} column
 'amplificationCall' has been appended to the segmentation table.
}

\section{The UCSF caller}{
  If \code{method == "ucsf-exp"}, then segments are called using [1], i.e.
  a segment is called an amplification if ...
}

\author{Henrik Bengtsson}

\references{
  [1] Fridlyand et al. \emph{Breast tumor copy number aberration
      phenotypes and genomic instability}, BMC Cancer, 2006. \cr
}

\seealso{
  \code{\link[PSCBS:callGainsAndLosses.CBS]{*callGainsAndLosses}()}.
  \code{\link[PSCBS:callOutliers.CBS]{*callOutliers}()}.
  For more information see \code{\link{CBS}}.
}


\keyword{internal}
\keyword{methods}
