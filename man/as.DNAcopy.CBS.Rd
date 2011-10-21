%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  DNAcopy.EXTS.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{as.DNAcopy.CBS}
\alias{as.DNAcopy.CBS}
\alias{CBS.as.DNAcopy}
\alias{as.DNAcopy.CBS}
\alias{as.DNAcopy,CBS-method}


\title{Coerces a CBS object to a DNAcopy object}

\description{
 Coerces a CBS object to a DNAcopy object.
}

\usage{\method{as.DNAcopy}{CBS}(fit, ...)}

\arguments{
  \item{fit}{A \code{\link{CBS}} object."}
  \item{...}{Not used.}
}

\value{
  Returns a \code{\link{DNAcopy}} object (of the \pkg{DNAcopy} package).
}

\examples{
 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Simulating copy-number data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
set.seed(0xBEEF)

# Number of loci
J <- 1000

mu <- double(J)
mu[200:300] <- mu[200:300] + 1
mu[350:400] <- NA # centromere
mu[650:800] <- mu[650:800] - 1
eps <- rnorm(J, sd=1/2)
y <- mu + eps
x <- sort(runif(length(y), max=length(y))) * 1e5
w <- runif(J)
w[650:800] <- 0.001


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Segmentation
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fit <- segmentByCBS(y, x=x)
print(fit)
plotTracks(fit)


 
# Coerce an CBS object to a DNAcopy object
fitD <- as.DNAcopy(fit)

# Coerce an DNAcopy object to a CBS object
fitC <- as.CBS(fitD)

# Sanity check
fitD2 <- as.DNAcopy(fit)
stopifnot(all.equal(fitD2, fitD))

fitC2 <- as.CBS(fitD2)
stopifnot(all.equal(fitC2, fitC))

}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  \code{\link[PSCBS:as.CBS.DNAcopy]{as.CBS()}}.
  For more information see \code{\link{CBS}}.
}


\keyword{internal}
\keyword{methods}
