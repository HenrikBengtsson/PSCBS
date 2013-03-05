%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  CBS.CALL.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{callOutliers.CBS}
\alias{callOutliers.CBS}
\alias{CBS.callOutliers}
\alias{callOutliers,CBS-method}

\title{Calls outliers}

\description{
 Calls outliers.
}

\usage{
\method{callOutliers}{CBS}(fit, adjust=1, method=c("ucsf-mad"), ...)
}

\arguments{
 \item{adjust}{A positive scale factor adjusting the sensitivity of the
   caller, where a value less (greater) than 1.0 makes the caller
   less (more) sensitive.}
 \item{method}{A \code{\link[base]{character}} string specifying the calling algorithm to use.}
 \item{...}{Additional/optional arguments used to override the default
   parameters used by the caller.}
}

\value{
 Returns a \code{\link[PSCBS]{CBS}} object where \code{\link[base]{logical}} columns
 'negOutlierCall' and 'posOutlierCall' have been appended
 to the segmentation table.
}

\section{The UCSF caller}{
  If \code{method == "ucsf-mad"}, then loci are called using [1];
 "Finally, to identify single technical or biological outliers such
  as high level amplifications, the presence of the outliers within
  a segment was allowed by assigning the original observed log2ratio
  to the clones for which the observed values were more than four
  tumor-specific MAD away from the smoothed values." [1; Suppl. Mat.]
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\references{
  [1] Fridlyand et al. \emph{Breast tumor copy number aberration
      phenotypes and genomic instability}, BMC Cancer, 2006. \cr
}

\seealso{
  \code{\link[PSCBS:callGainsAndLosses.CBS]{*callGainsAndLosses}()}.
  \code{\link[PSCBS:callAmplifications.CBS]{*callAmplifications}()}.
  For more information see \code{\link{CBS}}.
}


\keyword{internal}
\keyword{methods}
