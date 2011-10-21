%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  CBS.CALL.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 \name{callGainsAndLosses.CBS}
\alias{callGainsAndLosses.CBS}
\alias{CBS.callGainsAndLosses}
\alias{callGainsAndLosses.CBS}
\alias{callGainsAndLosses,CBS-method}


 \title{Calls gains and losses}

 \description{
  Calls gains and losses.
 }

 \usage{\method{callGainsAndLosses}{CBS}(fit, adjust=1, method=c("ucsf-mad"), ...)}

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
  'lossCall' and 'gainCall' have been appended to the segmentation table.
 }

 \section{The UCSF caller}{
   If \code{method == "ucsf-mad"}, then segments are called using [1], i.e.
   a segment is called gained or lost if its segment level is
   at least two standard deviations away from the median segment level
   on Chr1-22, where standard deviation is estimated using MAD.
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


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Segmentation with some known change points
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
knownSegments <- data.frame(
  chromosome=c(    0,   0),
  start     =x[c(  1, 401)],
  end       =x[c(349,   J)]+c(0,10e6)
)
fit2 <- segmentByCBS(y, x=x, knownSegments=knownSegments, verbose=TRUE)
print(fit2)
plotTracks(fit2)
abline(v=c(knownSegments$start, knownSegments$end)/1e6, lty=3)


# Segment just the (empty) centromere. All statistics will be NAs.
knownSegments <- data.frame(
  chromosome=c(    0),
  start     =x[c(350)],
  end       =x[c(400)]
)
fit3 <- segmentByCBS(y, x=x, knownSegments=knownSegments, verbose=TRUE)
print(fit3)



# If one specify the (empty) centromere as a segment, then its
# estimated statistics will be NAs.  [TODO /HB 2011-10-19]
knownSegments <- data.frame(
  chromosome=c(    0,   0,   0),
  start     =x[c(  1, 350, 401)],
  end       =x[c(349, 400,   J)]
)
fit4 <- segmentByCBS(y, x=x, knownSegments=knownSegments, verbose=TRUE)
print(fit4)
plotTracks(fit4)
abline(v=c(knownSegments$start, knownSegments$end)/1e6, lty=3)





  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# CALLS
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Call gains and losses by segments
fitC <- callGainsAndLosses(fit)

# Call amplifications by segments
fitC <- callAmplifications(fitC)

# Call outliers by loci
fitC <- callOutliers(fitC)

 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \references{
   [1] Fridlyand et al. \emph{Breast tumor copy number aberration
       phenotypes and genomic instability}, BMC Cancer, 2006. \cr
 }

 \seealso{
   \code{\link[PSCBS:callAmplifications.CBS]{*callAmplifications}()}.
   \code{\link[PSCBS:callOutliers.CBS]{*callOutliers}()}.
   For more information see \code{\link{CBS}}.
 }


\keyword{internal}
\keyword{methods}
