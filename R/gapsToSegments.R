###########################################################################/**
# @set "class=data.frame"
# @RdocMethod gapsToSegments
# @alias gapsToSegments
#
# @title "Gets the genomic segments that are complementary to the gaps"
#
# \description{
#  @get "title", with default chromosome boundaries being \code{-Inf}
#  and \code{+Inf}.
# }
#
# @synopsis
#
# \arguments{
#   \item{gaps}{A @data.frame with columns \code{chromosome}, \code{start},
#     and \code{stop}. Any overlapping gaps will throw an error.}
#   \item{resolution}{A non-negative @numeric specifying the minimum
#     length unit, which by default equals one nucleotide/base pair.}
#   \item{minLength}{Minimum length of segments to be kept.}
#   \item{dropGaps}{If @TRUE, the gaps themselves are not part of the output.}
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns @data.frame of least one row with columns \code{chromosome}
#   if that argument is given), \code{start}, \code{stop} and \code{length}.
#   The segments are ordered along the genome.
# }
#
# @author "HB"
#
# \seealso{
#   @see "findLargeGaps".
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("gapsToSegments", "data.frame", function(gaps, resolution=1L, minLength=0L, dropGaps=FALSE, ...) {
  # To please R CMD check
  chromosome <- NULL; rm(list="chromosome")

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'gaps':
  keys <- colnames(gaps)
  .stop_if_not(all(is.element(c("start", "end"), keys)))
  .stop_if_not(all(gaps$start <= gaps$end, na.rm=TRUE))
  hasChr <- is.element("chromosome", keys)

  ## Nothing more to do?
  if (nrow(gaps) == 0L) {
    knownSegments <- data.frame(chromosome=integer(1L), start=-Inf, end=+Inf)
    if (!hasChr) knownSegments$hromosome <- NULL
    return(knownSegments)
  }

  # Order gaps by the genome
  o <- order(gaps$chromosome, gaps$start, gaps$end)
  gaps <- gaps[o,]

  # For each chromosome...
  knownSegments <- NULL
  chromosomes <- sort(unique(gaps$chromosome))
  for (chr in chromosomes) {
    gapsCC <- subset(gaps, chromosome == chr)
    nCC <- nrow(gapsCC)

    starts <- gapsCC$start
    ends <- gapsCC$end

    # Assert that no overlapping gaps where specified
    if (!all(starts[-1] >= ends[-nCC], na.rm=TRUE)) {
      print(knownSegments)
      stop("INTERNAL ERROR: Detected overlapping gaps on chromosome ", chr, " in argument 'gaps'.")
    }

    # All boundaries in order
    # (this is possible because gaps are non-overlapping)
    naValue <- NA_real_
    if (dropGaps) {
      bps <- rep(naValue, times=2*nCC)
      bps[seq(from=1, to=2*nCC, by=2)] <- starts - resolution
      bps[seq(from=2, to=2*nCC, by=2)] <- ends + resolution
      bps <- c(-Inf, bps, +Inf)
      dim(bps) <- c(2L, nCC+1L)
    } else {
      bps <- rep(naValue, times=4*nCC)
      bps[seq(from=1, to=4*nCC, by=4)] <- starts - resolution
      bps[seq(from=2, to=4*nCC, by=4)] <- starts
      bps[seq(from=3, to=4*nCC, by=4)] <- ends
      bps[seq(from=4, to=4*nCC, by=4)] <- ends + resolution
      bps <- c(-Inf, bps, +Inf)
      dim(bps) <- c(2L, 2*nCC+1L)
    }

    knownSegmentsCC <- data.frame(chromosome=chr, start=bps[1L,], end=bps[2L,])

    knownSegments <- rbind(knownSegments, knownSegmentsCC)
  } # for (chr ...)

#  o <- with(knownSegments, order(chromosome, start, end))
#  knownSegments <- knownSegments[o,]
#  rownames(knownSegments) <- NULL

  # Append segment lengths
  knownSegments$length <- knownSegments$end - knownSegments$start

  # Drop too short segments
  keep <- (knownSegments$length >= minLength)
  knownSegments <- knownSegments[keep,]


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate generated 'knownSegments'
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  .stop_if_not(is.data.frame(knownSegments))
  .stop_if_not(nrow(knownSegments) >= 1L)
  for (chr in sort(unique(knownSegments$chromosome))) {
    dd <- subset(knownSegments, chromosome == chr)

    # Known segments must not overlap
    if (!all(dd$start[-1] >= dd$end[-nrow(dd)], na.rm=TRUE)) {
      stop("INTERNAL ERROR: Detected overlapping segments on chromosome ", chr, " in generated 'knownSegments'.")
    }
  }

  knownSegments
}) # gapsToSegments()
