###########################################################################/**
# @RdocFunction randomSeed
#
# @title "Sets and resets the .Random.seed in the global environment"
#
# \description{
#  @get "title".
# }
#
# \usage{
#  @usage randomSeed
# }
#
# \arguments{
#   \item{action}{A @character string specifying the action.}
#   \item{seed}{Random seed to be set; only for \code{action="set"}.
#     If \code{length(seed) == 1}, then \code{set.seed(seed)} is
#     used, otherwise \code{.Random.seed} is assigned the value.}
# }
#
# \value{
#   Returns invisibly the previous/current \code{.Random.seed}.
# }
#
# @author "HB"
#
# @keyword IO
#*/###########################################################################
randomSeed <- local({
  oldSeed <- NULL
  lecuyerSeed <- NULL

  genv <- globalenv()

  getSeed <- function() {
    if (exists(".Random.seed", envir=genv, inherits=FALSE)) {
      get(".Random.seed", envir=genv, inherits=FALSE)
    } else {
      NULL
    }
  }

  setSeed <- function(seed) {
    if (is.null(seed)) {
      if (exists(".Random.seed", envir=genv, inherits=FALSE))
        rm(list=".Random.seed", envir=genv, inherits=FALSE)
    } else {
      oldSeed <<- getSeed()
      if (length(seed) == 1L) {
        set.seed(seed)
      } else {
        assign(".Random.seed", seed, envir=genv, inherits=FALSE)
      }
    }
  }

  advanceSeed <- function() {
    ## Nothing to do?
    if (RNGkind()[1L] != "L'Ecuyer-CMRG") return()

    if (is.null(lecuyerSeed)) {
      stats::runif(1)
      lecuyerSeed <- getSeed()
    }

    lecuyerSeed <<- parallel::nextRNGStream(lecuyerSeed)
    setSeed(lecuyerSeed)
  }


  function(action=c("set", "advance", "reset", "get"), seed=NULL) {
    action <- match.arg(action)
    currSeed <- getSeed()

    if (action == "set") {
      setSeed(seed)
    } else if (action == "advance") {
      advanceSeed()
    } else if (action == "reset") {
      setSeed(oldSeed)
    } else if (action == "get") {
      return(currSeed)
    }

    invisible(currSeed)
  }
}) # randomSeed()
