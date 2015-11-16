rngSeeder <- local({
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


  function(action=c("set", "reset", "get"), seed=NULL) {
    action <- match.arg(action)
    currSeed <- getSeed()

    if (action == "set") {
      setSeed(seed)
      invisible(currSeed)
    } else if (action == "reset") {
      setSeed(oldSeed)
      invisible(currSeed)
    } else if (action == "get") {
      currSeed
    }
  }
}) # rngSeeder()
