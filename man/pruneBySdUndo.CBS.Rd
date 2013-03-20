%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  CBS.PRUNE.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{pruneBySdUndo.CBS}
\alias{pruneBySdUndo.CBS}
\alias{CBS.pruneBySdUndo}
\alias{pruneBySdUndo,CBS-method}

\title{Prune the CBS profile by dropping change points that are too small}

\description{
 Prune the CBS profile by dropping change points that are too small, where "too small" means that the amplitude of the
 change points is less than a multiple of the overall standard deviation
 of the copy-number signals.
}

\usage{
\method{pruneBySdUndo}{CBS}(fit, rho=3, sigma="DNAcopy", ..., verbose=FALSE)
}

\arguments{
  \item{fit}{A \code{\link{CBS}} object.}
  \item{rho}{A positive \code{\link[base]{double}} scalar specifying the number of standard
    deviations (\code{rho*sigma}) required in order to keep a change point.
    More change points are dropped the greater this value is.}
  \item{sigma}{The whole-genome standard deviation of the locus-level
    copy number signals.  The default is to calculate it from the data
    and as done in the \pkg{DNAcopy} package.}
  \item{...}{(Optional) Additional arguments passed to the standard
    deviation estimator function.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
  Returns a \code{\link{CBS}} object (of the same class as \code{fit}).
}

\details{
 This method corresponds to using the \code{undo} argument when calling
 \code{\link{segmentByCBS}}(), which in turn corresponds to using the
 \code{undo.splits="sdundo"} and \code{undo.SD} of the underlying
 \code{\link[DNAcopy]{segment}} method.
}

\examples{

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Simulating copy-number data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
set.seed(0xBEEF)

# Number of loci
J <- 1000

mu <- double(J)
mu[1:100] <- mu[1:100] + 0.3
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

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Post-segmentation pruning
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fitP <- pruneBySdUndo(fit, rho=1)
drawLevels(fitP, col="red")
}

\author{Henrik Bengtsson, Pierre Neuvial}


\keyword{internal}
\keyword{methods}
