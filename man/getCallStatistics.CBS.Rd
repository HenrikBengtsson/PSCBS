%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  CBS.CALL.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getCallStatistics.CBS}
\alias{getCallStatistics.CBS}
\alias{CBS.getCallStatistics}
\alias{getCallStatistics.CBS}
\alias{getCallStatistics,CBS-method}


\title{Calculates various call statistics per chromosome}

\description{
 Calculates various call statistics per chromosome.
}

\usage{\method{getCallStatistics}{CBS}(fit, ...)}

\arguments{
 \item{...}{Not used.}
}

\value{
 Returns a CxK \code{\link[base]{data.frame}}, where C is the number of chromosomes
 and (K-2)/2 is the number of call types.
 The first column is the chromosome index and the second is the
 length of the chromosome.
 The following columns contains call summaries per chromosome.
 For each chromosome and call type, the total length of such calls
 on that chromosome is reported together how large of a fraction
 of the chromosome such calls occupy.
}

\details{
  The estimators implemented here are based solely on the
  segmentation results, which is very fast.
  In the original proposal by Fridlyand et al. [1], the authors
  estimates the parameters by converting segment-level calls back
  to locus-level calls and there do the calculations.
  The difference between the two approaches should be minor,
  particularly for large density arrays.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\references{
  [1] Fridlyand et al. \emph{Breast tumor copy number aberration
      phenotypes and genomic instability}, BMC Cancer, 2006. \cr
}

\seealso{
  For more information see \code{\link{CBS}}.
}


\keyword{internal}
\keyword{methods}
