%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  PairedPSCBS.estimateDeltaAB.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{estimateDeltaABBySmallDH.PairedPSCBS}
\alias{estimateDeltaABBySmallDH.PairedPSCBS}
\alias{PairedPSCBS.estimateDeltaABBySmallDH}
\alias{estimateDeltaABBySmallDH.PairedPSCBS}
\alias{estimateDeltaABBySmallDH,PairedPSCBS-method}


\title{Estimate a threshold for calling allelic balance from DH}

\description{
 Estimate a threshold for calling allelic balance from DH.
}

\usage{\method{estimateDeltaABBySmallDH}{PairedPSCBS}(fit, q1=0.05, q2=0.9, ..., verbose=FALSE)}

\arguments{
  \item{q1}{A \code{\link[base]{numeric}} value specifying the weighted quantile of the
   segment-level DHs used to identify segments with small DH means.}
  \item{q2}{A \code{\link[base]{numeric}} value specifying the quantile of the locus-level
   DH signals for those segments with small DH mean levels.}
  \item{...}{Not used.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
  Returns the threshold estimate as a \code{\link[base]{numeric}} scalar.
}

\section{Algorithm}{
 \itemize{
  \item Grabs the segment-level DH estimates.
  \item Calculate segment weights proportional to the number
        of heterozygous SNPs.
  \item Calculate \eqn{\Delta} as the 5\% quantile of the weighted DH means.
  \item Choose the segments with means less than \eqn{\Delta}.
  \item Calculate threshold \eqn{\Delta_{AB}} as the 90\% "symmetric" quantile
        of the observed locus-level DHs from the selected segments
        in Step 4.
        The q:th "symmetric" quantile is estimated by estimating
        the ((1-q), 50\%) quantiles, calculating their distance as
        "50\%-(1-q)" and add to the median (50\%), i.e.
        "median + (median-(1-q))" = "2*median-1 + q", which should
        equal q if the distribution is symmetric.
 }
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  Instead of calling this method explicitly, it is recommended
  to use the \code{\link[PSCBS:estimateDeltaAB.PairedPSCBS]{*estimateDeltaAB}()} method.
}


\keyword{internal}
\keyword{methods}
