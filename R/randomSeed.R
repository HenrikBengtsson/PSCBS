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
#   \item{seed}{Random seed to be set; only for \code{action="set"}.}
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
  oldKind <- NULL

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
      set.seed(seed)
    }
  }


  function(action=c("set", "advance", "reset", "get"), seed=NULL) {
    action <- match.arg(action)
    currSeed <- getSeed()

    if (action == "set") {
      setSeed(seed)
      invisible(currSeed)
    } else if (action == "advance") {
      ## TO DO: Just a stub for now
      invisible(currSeed)
    } else if (action == "reset") {
      setSeed(oldSeed)
      invisible(currSeed)
    } else if (action == "get") {
      currSeed
    }
  }
}) # randomSeed()
