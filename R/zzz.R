.use <- function(name, package, ..., mode = "function") {
  get(name, mode = mode, envir = getNamespace(package))
}


.onLoad <- function(libname, pkgname) {
  ## covr: skip=5
  ns <- getNamespace(pkgname)
  pkg <- Package(pkgname)
  # Assign '.PSCBS' object [since 'PSCBS' is a constructor/Class].
  name <- sprintf(".%s", pkgname)
  assign(name, pkg, envir = ns, inherits = FALSE)
}


.onAttach <- function(libname, pkgname) {
  # Copy some pre-memoized CBS-parameter calculations to the 'R.cache'
  # cache.  This speeds up the calculation for common CBS use cases.
  .prememoize()

  # Get '.PSCBS' object [since 'PSCBS' is a constructor/Class].
  name <- sprintf(".%s", pkgname)
  pkg <- get(name, envir = getNamespace(pkgname), inherits = FALSE)
  startupMessage(pkg)
}
