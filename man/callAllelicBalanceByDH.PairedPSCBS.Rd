%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  PairedPSCBS.callAB.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{callAllelicBalanceByDH.PairedPSCBS}
\alias{callAllelicBalanceByDH.PairedPSCBS}
\alias{PairedPSCBS.callAllelicBalanceByDH}
\alias{callAllelicBalanceByDH,PairedPSCBS-method}

\title{Calls segments that are in allelic balance}

\description{
 Calls segments that are in allelic balance by thresholding on DH using a predetermined threshold.
 The variability of the DH mean levels is taken into account via a
 bootstrap estimator.
}

\usage{
\method{callAllelicBalanceByDH}{PairedPSCBS}(fit, delta=estimateDeltaAB(fit, flavor = "qq(DH)"), alpha=0.05, ..., verbose=FALSE)
}

\arguments{
  \item{flavor}{A \code{\link[base]{character}} string specifying which type of
   call to use.}
  \item{delta}{(Tuning parameter) A non-negative \code{\link[base]{numeric}} threshold.}
  \item{alpha}{A \code{\link[base]{numeric}} in [0,1] specifying the upper and lower
    quantiles calculated by the bootstrap estimator.}
  \item{...}{Additional arguments passed to the caller.}
}

\value{
  Returns a \code{\link{PairedPSCBS}} object with allelic-balance calls.
}

\author{Henrik Bengtsson}

\section{Algorithm}{
 \itemize{
   \item Foo
   \item Bar
 }
}

\seealso{
  Instead of calling this method explicitly, it is recommended
  to use the \code{\link[PSCBS:callAllelicBalance.PairedPSCBS]{*callAllelicBalance}()} method.
}


\keyword{internal}
\keyword{methods}
