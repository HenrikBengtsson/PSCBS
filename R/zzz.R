# Allows conflicts. For more information, see library() and
# conflicts() in [R] base.
.conflicts.OK <- TRUE

## .onAttach <- function(libname, pkgname) { # Only if NAMESPACE
.First.lib <- function(libname, pkgname) {
  pd <- utils::packageDescription(pkgname);

  # Inform user if DNAcopy is missing
  if (isPackageInstalled("DNAcopy")) {
    library("DNAcopy");
  } else {
    msg <- "The Bioconductor package 'DNAcopy' is not installed. Please see http://www.bioconductor.org/ on how to install it, or try calling installDNAcopy().";
    hrule <- paste(rep("*", times=getOption("width", 80L)-1L), collapse="");
    packageStartupMessage(sprintf("%s\nNOTE: %s\n%s", hrule, msg, hrule));
  }

  packageStartupMessage(pkgname, " v", pd$Version, " (", 
    pd$Date, ") successfully loaded. See ?", pkgname, " for help."); 
}
