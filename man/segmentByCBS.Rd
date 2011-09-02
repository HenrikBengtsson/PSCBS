%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  segmentByCBS.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{segmentByCBS}
\alias{segmentByCBS.default}
\alias{segmentByCBS}

 \title{Segment genomic signals using the CBS method}

 \description{
  Segment genomic signals using the CBS method of the \pkg{DNAcopy} package.
  This is a convenient low-level wrapper for the \code{DNAcopy::segment()}
  method.  It is intended to be applied to one sample and one chromosome
  at the time.
  For more details on the Circular Binary Segmentation (CBS) method
  see [1,2].
 }

 \usage{\method{segmentByCBS}{default}(y, chromosome=0, x=NULL, index=seq(along = y), w=NULL, undo=Inf, ..., joinSegments=TRUE, knownCPs=NULL, columnNamesFlavor=c("PSCBS", "DNAcopy"), seed=NULL, verbose=FALSE)}

 \arguments{
   \item{y}{A \code{\link[base]{numeric}} \code{\link[base]{vector}} of J genomic signals to be segmented.}
   \item{chromosome}{(Optional) An \code{\link[base]{integer}} scalar
       (or a \code{\link[base]{vector}} of length J contain a unique value).
       Only used for annotation purposes.}
   \item{x}{Optional \code{\link[base]{numeric}} \code{\link[base]{vector}} of J genomic locations.
            If \code{\link[base]{NULL}}, index locations \code{1:J} are used.}
   \item{w}{Optional \code{\link[base]{numeric}} \code{\link[base]{vector}} in [0,1] of J weights.}
   \item{undo}{A non-negative \code{\link[base]{numeric}}.  If less than +\code{\link[base:is.finite]{Inf}}, then
       arguments \code{undo.splits="sdundo"} and \code{undo.SD=undo}
       are passed to \code{DNAcopy::segment()}.}
   \item{...}{Additional arguments passed to the segmentation function.}
   \item{joinSegments}{If \code{\link[base:logical]{TRUE}}, there are no gaps between neighboring
     segments.
     If \code{\link[base:logical]{FALSE}}, the boundaries of a segment are defined by the support
     that the loci in the segments provides, i.e. there exist a locus
     at each end point of each segment.  This also means that there
     is a gap between any neighboring segments, unless the change point
     is in the middle of multiple loci with the same position.
     The latter is what \code{DNAcopy::segment()} returns.
   }
   \item{knownCPs}{Optional \code{\link[base]{numeric}} \code{\link[base]{vector}} of known
     change point locations.}
   \item{index}{An \code{\link[base]{integer}} \code{\link[base]{vector}} of length J specifying the
     genomewide indices of the loci.}
   \item{columnNamesFlavor}{A \code{\link[base]{character}} string specifying how column names
     of the returned data frame should be named.}
   \item{seed}{An (optional) \code{\link[base]{integer}} specifying the random seed to be
     set before calling the segmentation method.  The random seed is
     set to its original state when exiting.  If \code{\link[base]{NULL}}, it is not set.}
   \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
 }

 \value{
   Returns the fit object.
 }

 \details{
   Internally \code{\link[DNAcopy]{segment}} is used to segment the signals.
   This segmentation method support weighted segmentation.

   The "DNAcopy::segment" implementation of CBS uses approximation
   through random sampling for some estimates.  Because of this,
   repeated calls using the same signals may result in slightly
   different results, unless the random seed is set/fixed.
 }

 \section{Missing and non-finite values}{
   Signals as well as genomic positions may contain missing
   values, i.e. \code{\link[base]{NA}}s or \code{\link[base:is.finite]{NaN}}s.  If so, they are silently excluded
   before performing the segmentation.

   None of the input signals may have infinite values,
   i.e. -\code{\link[base:is.finite]{Inf}} or +\code{\link[base:is.finite]{Inf}}. If so, an informative error is thrown.
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
mu[650:800] <- mu[650:800] - 1
eps <- rnorm(J, sd=1/2)
y <- mu + eps
x <- sort(runif(length(y), max=length(y))) * 1e5
w <- runif(J)
w[650:800] <- 0.001


xlab <- "Position (Mb)"
ylim <- c(-3,3)
xMb <- x/1e6
plot(xMb,y, pch=20, col="#aaaaaa", xlab=xlab, ylim=ylim)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Segment
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fit <- segmentByCBS(y, x=x)
print(fit)
drawLevels(fit, col="red", lwd=2, xScale=1e-6)



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# TESTS
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fit <- segmentByCBS(y, x=x, seed=0xBEEF)
print(fit)
##   id chromosome       start      end nbrOfLoci    mean
## 1  y          0    55167.82 20774251       201  0.0164
## 2  y          0 20774250.85 29320105        99  1.0474
## 3  y          0 29320104.86 65874675       349 -0.0227
## 4  y          0 65874675.06 81348129       151 -1.0813
## 5  y          0 81348129.20 99910827       200 -0.0612


# Test #1: Reverse the ordering and segment
fitR <- segmentByCBS(rev(y), x=rev(x), seed=0xBEEF)
# Sanity check
stopifnot(all.equal(fitR$output, fit$output))

# Test #2: Reverse, but preserve ordering of 'data' object
fitRP <- segmentByCBS(rev(y), x=rev(x), preserveOrder=TRUE)
stopifnot(all.equal(fitRP$output, fit$output))


# (Test #3: Change points inbetween data points at the same locus)
x[650:654] <- x[649]
fitC <- segmentByCBS(rev(y), x=rev(x), preserveOrder=TRUE, seed=0xBEEF)

# Test #4: Allow for some missing values in signals
y[450] <- NA
fitD <- segmentByCBS(y, x=x, seed=0xBEEF)


# Test #5: Allow for some missing genomic annotations
x[495] <- NA
fitD <- segmentByCBS(y, x=x, seed=0xBEEF)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# MISC.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Emulate a centromere
x[650:699] <- NA
fit <- segmentByCBS(y, x=x, seed=0xBEEF)
xMb <- x/1e6
plot(xMb,y, pch=20, col="#aaaaaa", xlab=xlab, ylim=ylim)
drawLevels(fit, col="red", lwd=2, xScale=1e-6)

fitC <- segmentByCBS(y, x=x, joinSegments=FALSE, seed=0xBEEF)
drawLevels(fitC, col="blue", lwd=2, xScale=1e-6)



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Multiple chromosomes
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Appending CBS results
fit1 <- segmentByCBS(y, chromosome=1, x=x)
fit2 <- segmentByCBS(y, chromosome=2, x=x)
fit <- append(fit1, fit2)
print(fit)
plotTracks(fit, subset=NULL, pch=20, lwd=2, Clim=c(-3,3))


# Segmenting multiple chromosomes at once
chromosomeWG <- rep(1:2, each=J)
xWG <- rep(x, times=2)
yWG <- rep(y, times=2)
fitWG <- segmentByCBS(yWG, chromosome=chromosomeWG, x=xWG)
print(fitWG)
plotTracks(fitWG, subset=NULL, pch=20, lwd=2, Clim=c(-3,3))

# Assert same results
fit$data[,"index"] <- fitWG$data[,"index"] # Ignore 'index'
stopifnot(all.equal(fitWG$data, fit$data))
stopifnot(all.equal(fitWG$output, fit$output))
}

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \references{
  [1] A.B. Olshen, E.S. Venkatraman (aka Venkatraman E. Seshan),
      R. Lucito and M. Wigler, \emph{Circular binary segmentation for
      the analysis of array-based DNA copy number data},
      Biostatistics, 2004.\cr
  [2] E.S. Venkatraman and A.B. Olshen, \emph{A faster circular binary
      segmentation algorithm for the analysis of array CGH data}.
      Bioinformatics, 2007.\cr
 }


\keyword{IO}
