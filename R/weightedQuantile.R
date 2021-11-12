############################################################################/**
# @RdocDefault weightedQuantile
#
# @title "Weighted Quantile Value"
#
# @synopsis
#
# \description{
#   Computes a weighted quantile of a numeric vector.
# }
#
# \arguments{
#   \item{x}{a @numeric @vector containing the values whose weighted
#            quantile is to be computed.}
#   \item{w}{a numeric @vector of weights the same length as
#            \code{x} giving the weights to use for each element of \code{x}.
#            Negative weights are treated as zero weights.
#            Default value is equal weight to all values.}
#   \item{probs}{a @numeric @vector of quantiles in [0,1] to be retrieved.}
#   \item{na.rm}{a @logical value indicating whether @NA values in
#            \code{x} should be stripped before the computation proceeds,
#            or not.}
#   \item{method}{If \code{"wtd.quantile"}, then an internal copy of
#            \code{Hmisc::wtd.quantile()} is used.
#            No other methods are currently supported.}
#   \item{...}{Additional arguments passed to the estimator.}
# }
#
# \value{
#   Returns the weighted quantile.
# }
#
# @author "HB"
#
# \seealso{
#   Internally the following functions may be used:
#   @see "stats::quantile" (if no weights are specified), or an internal
#   copy of \code{Hmisc::wtd.quantile()}.
#   For a weighted median estimator, @see "matrixStats::weightedMedian"
#   of the \pkg{matrixStats} package.
# }
#
# @keyword univar
# @keyword robust
# @keyword internal
#*/############################################################################
setMethodS3("weightedQuantile", "default", function(x, w, probs=c(0, 0.25, 0.5, 0.75, 1), na.rm=TRUE, method=c("wtd.quantile"), ...) {
  # Argument 'x':
  x <- Arguments$getNumerics(x)

  # Argument 'w':
  if (missing(w)) {
    # By default use weights that are one.
    w <- rep(1, times=length(x))
  } else {
    w <- Arguments$getNumerics(w, range=c(0,Inf), length=rep(length(x), times=2L))
  }

  naValue <- NA
  storage.mode(naValue) <- storage.mode(x)

  # Argument 'na.rm':
  if (is.na(na.rm)) {
    # There are no NAs
  } else if (isTRUE(na.rm)) {
    # Remove values that are NA's
    tmp <- !(is.na(x) | is.na(w))
    x <- .subset(x, tmp)
    w <- .subset(w, tmp)
  } else if (anyNA(x)) {
    return(naValue)
  }

  # Argument 'method':
  method <- match.arg(method)

  # Remove values with zero (and negative) weight. This will:
  # (1) take care of the case when all weights are zero,
  # (2) it will most likely speed up the sorting.
  n <- length(w)
  tmp <- (w > 0)
  if (!all(tmp)) {
    x <- .subset(x, tmp)
    w <- .subset(w, tmp)
    n <- sum(tmp)
  }

  # Are there any values left to calculate the weighted median of?
  if (n == 0) {
    return(naValue)
  } else if (n == 1) {
    return(x)
  }

  # Are any weights Inf? Then treat them with equal weight and all others
  # with weight zero. If they have equal weight, regular quantile
  # can be used instead, which is assumed to be faster.
  tmp <- is.infinite(w)
  if (any(tmp)) {
    x <- .subset(x, tmp)
    # Here we know there are no NAs.
    return(quantile(x, probs=probs, na.rm=FALSE, ...))
  }

  # Here we know that there are no missing values in the data
  .wtd.quantile(x, weights=w, probs=probs)
}) # weightedQuantile()



## The wtd.quantile() function originates from Hmisc 4.6-0 (2021-10-07), which
## was released under GPL (>= 2).  The reasons for copying (and pruning) it
## instead of adding 'Hmisc' as a dependency are several: (i) Hmisc has a large
## number of dependencies, (ii) Hmisc requires R (>= 3.6.0) whereas we try to
## stay backward compatible with R (>= 3.2.0), and, most importantly, (iii) due
## to it's many dependencies Hmisc does not install out of the box on all
## systems, e.g. as of 2021-10-22 it is not available for macOS running on the
## M1 chip.
##
## CHANGES MADE:
## * wtd.quantile(): Dropped argument 'type' - always type="quantile"
## * wtd.quantile(): Dropped argument 'normwt' - always normwt=TRUE
## * wtd.quantile(): Dropped argument 'na.rm' - always na.rm=FALSE
.wtd.quantile <- function(x, weights, probs=c(0, .25, .5, .75, 1)) {
  if(any(probs < 0 | probs > 1))
    stop("Probabilities must be between 0 and 1 inclusive")

  ## NOTE: Data points with zero or NA weights have already been dropped
  ## by weightedQuantile() before calling this function

  ## Normalize weights
  weights <- weights/sum(weights)

  ## Order (x,weights) by x
  i <- order(x)
  x <- x[i]
  weights <- weights[i]

  nx <- length(x)

  ## Merge replicated 'x':s into single ones by combining their weights
  if (anyDuplicated(x)) {
    weights <- tapply(weights, INDEX = x, FUN = sum)
    ## The names of 'weights' holds the unique 'x' values
    xs <- names(weights)
    names(weights) <- NULL
    if (length(xs) == 0) stop("program logic error")
    x <- as.numeric(xs)
  }
  
  weights <- nx * weights
  cweights <- cumsum(weights)
  
  n     <- cweights[length(weights)]
  order <- 1 + (n - 1) * probs
  low   <- pmax(floor(order), 1)
  high  <- pmin(low + 1, n)
  order <- order %% 1
  ## Find low and high order statistics
  ## These are minimum values of x such that the cum. freqs >= c(low,high)
  allq <- approx(cweights, x, xout=c(low,high), 
                 method="constant", f=1, rule=2L)$y
  k <- length(probs)
  quantiles <- (1 - order)*allq[1:k] + order*allq[-(1:k)]

  ## Add 'probs' names
  digits <- if (k > 1) 2 - log10(diff(range(probs))) else 2
  names <- paste(format(round(100 * probs, digits = digits)), "%", sep = "")
  names(quantiles) <- names
  
  quantiles
}
