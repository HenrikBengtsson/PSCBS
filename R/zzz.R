## .onAttach <- function(libname, pkgname) { # Only if NAMESPACE
.First.lib <- function(libname, pkgname) {
  pd <- utils::packageDescription(pkgname);
  packageStartupMessage(pkgname, " v", pd$Version, " (", 
    pd$Date, ") successfully loaded. See ?", pkgname, " for help."); 
}
