# Allows conflicts. For more information, see library() and
# conflicts() in [R] base.
.conflicts.OK <- TRUE


.onAttach <- function(libname, pkgname) {
  pkg <- Package(pkgname);

  # Inform user if DNAcopy is missing
  if (isPackageInstalled("DNAcopy")) {
    # To please/trick R CMD check
    loadPackage <- base::library;
    loadPackage("DNAcopy");
  } else {
    msg <- "The Bioconductor package 'DNAcopy' is not installed. Please see http://www.bioconductor.org/ on how to install it, or try calling installDNAcopy().";
    hrule <- paste(rep("*", times=getOption("width", 80L)-1L), collapse="");
    packageStartupMessage(sprintf("%s\nNOTE: %s\n%s", hrule, msg, hrule));
  }

  startupMessage(pkg);
}


############################################################################
# HISTORY:
# 2011-07-23
# o Added a namespace to the package.
############################################################################
