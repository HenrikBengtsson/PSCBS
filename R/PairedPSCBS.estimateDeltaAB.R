###########################################################################/**
# @set class=PairedPSCBS
# @RdocMethod estimateDeltaAB
#
# @title "Estimate a threshold for calling allelic balance from DH"
#
# \description{
#  @get "title" to be used by the @seemethod "callAB" method.
# }
#
# @synopsis
#
# \arguments{
#   \item{scale}{An optional @numeric scale factor.}
#   \item{flavor}{A @character string specifying which type of
#    estimator to use.}
#   \item{...}{Additional arguments passed to the estimator.}
#   \item{max}{(Optional) The maximum estimate allowed. If greater than 
#    this value, the estimate will be truncated.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   Returns the threshold estimate as a @numeric scalar.
# }
#
# @author "HB"
#
# \seealso{
#   Internally, one of the following methods are used:
#   @seemethod "estimateDeltaABBySmallDH",
#   @seemethod "estimateStdDevForHeterozygousBAF",
#   @seemethod "estimateMeanForDH", and
#   @seemethod "estimateHighDHQuantileAtAB".
# }
#
#*/###########################################################################  
setMethodS3("estimateDeltaAB", "PairedPSCBS", function(this, scale=NULL, flavor=c("qq(DH)", "q(DH)", "mad(hBAF)", "median(DH)"), ..., max=Inf, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'flavor':
  flavor <- match.arg(flavor)

  # Argument 'max':
  max <- Arguments$getDouble(max, range=c(0,Inf))

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    pushState(verbose)
    on.exit(popState(verbose))
  }

  verbose && enter(verbose, "Estimating DH threshold for calling allelic imbalances")
  verbose && cat(verbose, "flavor: ", flavor)

  if (flavor == "mad(hBAF)") {
    if (is.null(scale)) scale <- 3
    verbose && cat(verbose, "scale: ", scale)
    # sigma = mad(hBAF) = 1.4826*median(|hBAF-m|),
    # where m = median(hBAF) ~= 1/2
    sd <- estimateStdDevForHeterozygousBAF(this, ..., verbose=verbose)
    verbose && printf(verbose, "sd: %.3g\n", sd)
    delta <- scale * sd
  } else if (flavor == "median(DH)") {
    if (is.null(scale)) scale <- 3
    verbose && cat(verbose, "scale: ", scale)
    # sigma = 1/2*1.4826*median(|hBAF-1/2|), 
    # because DH = 2*|hBAF-1/2|
    mu <- estimateMeanForDH(this, ..., verbose=verbose)
    verbose && printf(verbose, "mu: %.3g\n", mu)
    sd <- 1/2 * 1.4826 * mu
    verbose && printf(verbose, "sd: %.3g\n", sd)
    delta <- scale * sd
  } else if (flavor == "q(DH)") {
    if (is.null(scale)) scale <- 1
    verbose && cat(verbose, "scale: ", scale)
    delta <- estimateHighDHQuantileAtAB(this, scale=scale, ..., verbose=verbose)
  } else if (flavor == "qq(DH)") {
    if (is.null(scale)) scale <- 1
    verbose && cat(verbose, "scale: ", scale)
    delta <- estimateDeltaABBySmallDH(this, ..., verbose=verbose)
    delta <- scale * delta
  } else {
    stop("Unkown flavor: ", flavor)
  }

##   } else if (flavor == "DHskew") {
##     fit <- this
##     if (is.null(fit$output$dhSkew)) {
##       verbose && enter(verbose, "Estimating DH skewness for each segment")
##       fit <- applyByRegion(fit, FUN=.addTcnDhStatitics, verbose=less(verbose, 5))
##       verbose && exit(verbose)
##     }
##     mu <- fit$output$dhMean
##     skew <- fit$output$dhSkew
## 
##     deltaSkew <- -0.55
##     keep <- which(skew < deltaSkew)
##     verbose && printf(verbose, "Number of segments heavily skewed (< %.3f): %d\n", deltaSkew, length(keep))
##     # Sanity check
##     if (length(keep) == 0) {
##       stop("Cannot estimate DH threshold for AB. No segments with strong skewness exists.")
##     }
##     deltaDH <- median(mu[keep], na.rm=TRUE)
##     verbose && printf(verbose, "deltaDH: %.3g\n", deltaDH)
##     deltaDH <- 1.10*deltaDH
##     verbose && printf(verbose, "Adjusted +10%% deltaDH: %.3g\n", deltaDH)
## 
##     # sigma = 1/2*1.4826*median(|hBAF-1/2|), 
##     # because DH = 2*|hBAF-1/2|
##     mu <- estimateMeanForDH(this, delta=deltaDH, ...)
##     verbose && printf(verbose, "mu: %.3g\n", mu)
##     sd <- 1/2 * 1.4826 * mu
##     verbose && printf(verbose, "sd: %.3g\n", sd)
##  }

  verbose && printf(verbose, "Estimated delta: %.3g\n", delta)

  # Truncate estimate?
  if (delta > max) {
    warning("Estimated delta (%.3g) was greater than the maximum allowed value (%.3g).  The latter will be used instead.", delta, max)
    delta <- max
    verbose && printf(verbose, "Max delta: %.3g\n", max)
    verbose && printf(verbose, "Truncated delta: %.3g\n", delta)
  }

  verbose && exit(verbose)

  delta
}) # estimateDeltaAB()



setMethodS3("estimateStdDevForHeterozygousBAF", "PairedPSCBS", function(this, deltaDH=0.20, deltaTCN=5, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'deltaDH':
  deltaDH <- Arguments$getDouble(deltaDH, range=c(0,1))

  # Argument 'deltaTCN':
  deltaTCN <- Arguments$getDouble(deltaTCN, range=c(0,Inf))

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    pushState(verbose)
    on.exit(popState(verbose))
  }


  verbose && enter(verbose, "Estimating standard deviation of tumor BAFs for heterozygous SNPs")
  verbose && cat(verbose, "DH threshold: ", deltaDH)
  verbose && cat(verbose, "TCN threshold: ", deltaTCN)

  segs <- as.data.frame(this)

  verbose && cat(verbose, "Number of segments: ", nrow(segs))


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Find segments to be used for the estimation
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Find segments that have low DHs
  idxsDH <- which(segs$dhMean <= deltaDH)
  verbose && cat(verbose, "Identified segments with small DH levels: ", length(idxsDH))
  verbose && str(verbose, idxsDH)

  # Sanity check
  if (length(idxsDH) == 0) {
    stop("Cannot estimate standard deviation.  There exist no segments with DH less or equal to the given threshold: ", deltaDH)
  }

  # Find segments that have low TCNs
  idxsTCN <- which(segs$tcnMean <= deltaTCN)
  verbose && cat(verbose, "Identified segments with small TCN levels: ", length(idxsTCN))
  verbose && str(verbose, idxsTCN)

  # Sanity check
  if (length(idxsTCN) == 0) {
    stop("Cannot estimate standard deviation.  There exist no segments with TCN less or equal to the given threshold: ", deltaTCN)
  }

  # Segments with small DH and small TCN
  idxs <- intersect(idxsDH, idxsTCN)
  verbose && cat(verbose, "Identified segments with small DH and small TCN levels: ", length(idxs))
  verbose && str(verbose, idxs)

  # Sanity check
  if (length(idxs) == 0) {
    stop("Cannot estimate standard deviation.  There exist no segments with small DH and small TCN.")
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract data and estimate parameters
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract those segments
  verbose && enter(verbose, "Extracting identified segments")
  fitT <- extractRegions(this, idxs)
  verbose && exit(verbose)

  # Get the tumor BAFs for the heterozygous SNPs
  verbose && enter(verbose, "Extracting BAFs for the heterozygous SNPs")
  beta <- with(fitT$data, betaTN[muN == 1/2])
  verbose && str(verbose, beta)
  verbose && exit(verbose)

  # Estimate the standard deviation for those
  sd <- mad(beta, na.rm=TRUE)
  verbose && cat(verbose, "Estimated standard deviation: ", sd)


  verbose && exit(verbose)

  sd
}, private=TRUE) # estimateStdDevForHeterozygousBAF()




setMethodS3("estimateMeanForDH", "PairedPSCBS", function(this, deltaDH=0.20, deltaTCN=5, robust=TRUE, trim=0.05, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'deltaDH':
  deltaDH <- Arguments$getDouble(deltaDH, range=c(0,1))

  # Argument 'deltaTCN':
  deltaTCN <- Arguments$getDouble(deltaTCN, range=c(0,Inf))

  # Argument 'robust':
  robust <- Arguments$getLogical(robust)

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    pushState(verbose)
    on.exit(popState(verbose))
  }


  verbose && enter(verbose, "Estimating mean of tumor DHs for heterozygous SNPs")
  verbose && cat(verbose, "DH threshold: ", deltaDH)
  verbose && cat(verbose, "TCN threshold: ", deltaTCN)
  verbose && cat(verbose, "Robust estimator: ", robust)
  verbose && cat(verbose, "Trim: ", trim)

  segs <- as.data.frame(this)

  verbose && cat(verbose, "Number of segments: ", nrow(segs))

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Find segments to be used for the estimation
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Find segments that have low DHs
  idxsDH <- which(segs$dhMean <= deltaDH)
  verbose && cat(verbose, "Identified segments with small DH levels: ", length(idxsDH))
  verbose && str(verbose, idxsDH)

  # Sanity check
  if (length(idxsDH) == 0) {
    stop("Cannot estimate standard deviation.  There exist no segments with DH less or equal to the given threshold: ", deltaDH)
  }

  # Find segments that have low TCNs
  idxsTCN <- which(segs$tcnMean <= deltaTCN)
  verbose && cat(verbose, "Identified segments with small TCN levels: ", length(idxsTCN))
  verbose && str(verbose, idxsTCN)

  # Sanity check
  if (length(idxsTCN) == 0) {
    stop("Cannot estimate standard deviation.  There exist no segments with TCN less or equal to the given threshold: ", deltaTCN)
  }

  # Segments with small DH and small TCN
  idxs <- intersect(idxsDH, idxsTCN)
  verbose && cat(verbose, "Identified segments with small DH and small TCN levels: ", length(idxs))
  verbose && str(verbose, idxs)

  # Sanity check
  if (length(idxs) == 0) {
    stop("Cannot estimate standard deviation.  There exist no segments with small DH and small TCN.")
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract data and estimate parameters
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract those segments
  verbose && enter(verbose, "Extracting identified segments")
  fitT <- extractRegions(this, idxs)
  verbose && exit(verbose)

  # Get the tumor DHs for the heterozygous SNPs
  verbose && enter(verbose, "Extracting DHs for the heterozygous SNPs")
  rho <- with(fitT$data, rho[muN == 1/2])
  verbose && str(verbose, rho)
  verbose && exit(verbose)

  # Estimate the average for those
  rho <- rho[is.finite(rho)]
  if (robust) {
    mu <- median(rho, na.rm=FALSE)
    qlow <- quantile(rho, probs=0.05, na.rm=FALSE)
    delta <- mu-qlow
    print(list(qlow=qlow, mu=mu, delta=delta, "mu+delta"=mu+delta))
  } else {
    mu <- mean(rho, trim=trim, na.rm=FALSE)
  }
  verbose && cat(verbose, "Estimated mean: ", mu)


  verbose && exit(verbose)

  mu
}, private=TRUE) # estimateMeanForDH()



setMethodS3("estimateHighDHQuantileAtAB", "PairedPSCBS", function(this, quantile=0.99, scale=1, deltaDH=0.20, deltaTCN=5, robust=TRUE, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'quantile':
  quantile <- Arguments$getDouble(quantile, range=c(0.5,1))

  # Argument 'scale':
  scale <- Arguments$getDouble(scale, range=c(0,Inf))

  # Argument 'deltaDH':
  deltaDH <- Arguments$getDouble(deltaDH, range=c(0,1))

  # Argument 'deltaTCN':
  deltaTCN <- Arguments$getDouble(deltaTCN, range=c(0,Inf))

  # Argument 'robust':
  robust <- Arguments$getLogical(robust)

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    pushState(verbose)
    on.exit(popState(verbose))
  }


  verbose && enter(verbose, "Estimating DH quantile of tumor DHs for heterozygous SNPs")
  verbose && cat(verbose, "DH threshold: ", deltaDH)
  verbose && cat(verbose, "TCN threshold: ", deltaTCN)
  verbose && cat(verbose, "Robust estimator: ", robust)
  verbose && cat(verbose, "Scale factor: ", scale)

  segs <- as.data.frame(this)

  verbose && cat(verbose, "Number of segments: ", nrow(segs))


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Find segments to be used for the estimation
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Finding some segments that are likely to in allelic balance (AB)")

  # Find some segments that have low DHs
  idxsDH <- which(segs$dhMean <= deltaDH)
  verbose && cat(verbose, "Identified segments with small DH levels: ", length(idxsDH))
  verbose && str(verbose, idxsDH)

  # Sanity check
  if (length(idxsDH) == 0) {
    stop("Cannot estimate standard deviation.  There exist no segments with DH less or equal to the given threshold: ", deltaDH)
  }

  # Find segments that have low TCNs
  idxsTCN <- which(segs$tcnMean <= deltaTCN)
  verbose && cat(verbose, "Identified segments with small TCN levels: ", length(idxsTCN))
  verbose && str(verbose, idxsTCN)

  # Sanity check
  if (length(idxsTCN) == 0) {
    stop("Cannot estimate standard deviation.  There exist no segments with TCN less or equal to the given threshold: ", deltaTCN)
  }

  # Segments with small DH and small TCN
  idxs <- intersect(idxsDH, idxsTCN)
  verbose && cat(verbose, "Identified segments with small DH and small TCN levels: ", length(idxs))
  verbose && str(verbose, idxs)

  # Sanity check
  if (length(idxs) == 0) {
    stop("Cannot estimate standard deviation.  There exist no segments with small DH and small TCN.")
  }

  # Extract those segments
  verbose && enter(verbose, "Extracting identified segments")
  fitT <- extractRegions(this, idxs)
  verbose && exit(verbose)

  verbose && exit(verbose)



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract data and estimate parameters
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get the tumor DHs for the heterozygous SNPs
  verbose && enter(verbose, "Extracting DHs for the heterozygous SNPs")
  rho <- with(fitT$data, rho[muN == 1/2])
  verbose && str(verbose, rho)
  verbose && exit(verbose)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Estimating the DH quantile
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Estimating the quantile of interest")
  verbose && cat(verbose, "Quantile: ", quantile)

  # Drop missing values
  rho <- rho[is.finite(rho)]

  if (robust) {
    lq <- quantile(rho, probs=1-quantile, na.rm=FALSE)
    verbose && printf(verbose, "Estimated lower quantile (%.3f): %f\n", 1-quantile, lq)
    mu <- median(rho, na.rm=FALSE)
    verbose && cat(verbose, "Estimated median: ", mu)
    delta <- mu-lq
    verbose && printf(verbose, "Estimated \"spread\": %f\n", delta)
    uq <- mu + scale*delta
    verbose && printf(verbose, "Scale parameter: %f\n", scale)
    qs <- c(lq, mu, mu+delta, uq)
    names(qs) <- sprintf("%.1f%%", 100*c(1-quantile, 0.5, quantile, 0.5+scale*(quantile-0.5)))
    names(qs)[3:4] <- sprintf("%s*", names(qs)[3:4])
    attr(uq, "quantiles") <- qs
  } else {
    uq <- quantile(rho, probs=quantile, na.rm=FALSE)
  }

  names(uq) <- uq
  verbose && printf(verbose, "Estimated upper quantile (%.3f): %f\n", quantile, uq)

  verbose && exit(verbose)

  verbose && exit(verbose)

  uq
}, private=TRUE) # estimateHighDHQuantileAtAB()




###########################################################################/**
# @set class=PairedPSCBS
# @RdocMethod estimateDeltaABBySmallDH
#
# @title "Estimate a threshold for calling allelic balance from DH"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{q1}{A @numeric value specifying the weighted quantile of the
#    segment-level DHs used to identify segments with small DH means.}
#   \item{q2}{A @numeric value specifying the quantile of the locus-level
#    DH signals for those segments with small DH mean levels.}
#   \item{...}{Not used.}
#   \item{verbose}{See @see "R.utils::Verbose".} 
# }
#
# \value{
#   Returns the threshold estimate as a @numeric scalar.
# }
#
# \section{Algorithm}{
#  \itemize{
#   \item Grabs the segment-level DH estimates.
#   \item Calculate segment weights proportional to the number 
#         of heterozygous SNPs.
#   \item Calculate \eqn{\Delta} as the 5\% quantile of the weighted DH means.
#   \item Choose the segments with means less than \eqn{\Delta}.
#   \item Calculate threshold \eqn{\Delta_{AB}} as the 90\% "symmetric" quantile 
#         of the observed locus-level DHs from the selected segments 
#         in Step 4.
#         The q:th "symmetric" quantile is estimated by estimating 
#         the ((1-q), 50\%) quantiles, calculating their distance as
#         "50\%-(1-q)" and add to the median (50\%), i.e.
#         "median + (median-(1-q))" = "2*median-1 + q", which should
#         equal q if the distribution is symmetric.
#  }
# }
#
# @author "HB"
#
# \seealso{
#   Instead of calling this method explicitly, it is recommended
#   to use the @seemethod "estimateDeltaAB" method.
# }
#
# @keyword internal
#*/###########################################################################
setMethodS3("estimateDeltaABBySmallDH", "PairedPSCBS", function(fit, q1=0.05, q2=0.90, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
  # Argument 'q1' & 'q2':
  q1 <- Arguments$getDouble(q1, range=c(0,1))
  q2 <- Arguments$getDouble(q2, range=c(0,1))

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    pushState(verbose)
    on.exit(popState(verbose))
  }

  verbose && enter(verbose, "Estimating DH threshold for AB caller")
  verbose && cat(verbose, "quantile #1: ", q1)
  verbose && cat(verbose, "Symmetric quantile #2: ", q2)

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Extract the region-level estimates
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  segs <- getSegments(fit)
  dh <- segs$dhMean
  .stop_if_not(!is.null(dh))
  n <- segs$dhNbrOfLoci

  # Drop missing values
  keep <- (!is.na(dh) & !is.na(n))
  idxs <- which(keep)
  dh <- dh[idxs]
  n <- n[idxs]
  verbose && cat(verbose, "Number of segments: ", length(idxs))
  # Sanity check
  .stop_if_not(length(idxs) > 0)

  # Calculated weighted quantile
  weights <- n / sum(n)
  deltaDH <- weightedQuantile(dh, w=weights, probs=q1)
  verbose && printf(verbose, "Weighted %g%% quantile of DH: %f\n", 100*q1, deltaDH)

  # Identify segments with DH this small
  keep <- (dh <= deltaDH)
  idxs <- idxs[keep]
  verbose && cat(verbose, "Number of segments with small DH: ", length(idxs))
  # Sanity check
  .stop_if_not(length(idxs) > 0)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Extract the locus-level estimates
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Extract the regions of interest
  fitT <- extractRegions(fit, idxs)

  # Extract the data
  data <- fitT$data
  rho <- data$rho
  .stop_if_not(!is.null(rho))

  verbose && cat(verbose, "Number of data points: ", length(rho))

  # Drop missing values
  rho <- rho[is.finite(rho)]
  verbose && cat(verbose, "Number of finite data points: ", length(rho))

  qs <- quantile(rho, probs=c(1-q2, 1/2), na.rm=FALSE, names=FALSE)
  verbose && printf(verbose, "Estimate of (1-%.3g):th and 50%% quantiles: (%g,%g)\n", q2, qs[1], qs[2])
  deltaAB <- qs[2] + (qs[2]-qs[1])
  verbose && printf(verbose, "Estimate of %.3g:th \"symmetric\" quantile: %g\n", q2, deltaAB)
  

  # Sanity check
  deltaAB <- Arguments$getDouble(deltaAB)

  deltaAB
}, protected=TRUE) # estimateDeltaABBySmallDH()
