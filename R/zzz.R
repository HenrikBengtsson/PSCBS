# Allows conflicts. For more information, see library() and
# conflicts() in [R] base.
.conflicts.OK <- TRUE


.requirePkg <- function(name, quietly=FALSE) {
  # Nothing to do?
  if (is.element(sprintf("package:%s", name), search())) {
    return(invisible(TRUE));
  }
  if (quietly) {
    # Load the package (super quietly)
    res <- suppressPackageStartupMessages(require(name, character.only=TRUE, quietly=TRUE));
  } else {
    res <- require(name, character.only=TRUE);
  }
  if (!res) throw("Package not loaded: ", name);
  invisible(res);
} # .requirePkg()


# This function exists solely for the purpose of backward compatibility
# /HB 2013-09-27
.useAromaLight <- function(name=NULL, mode="function") {
  if (!isPackageInstalled("aroma.light")) {
    throw("Bioconductor package not installed: aroma.light");
  }

  # Are S3 method declared in NAMESPACE?
  ver <- packageVersion("aroma.light");
  properS3methods <- (ver >= "1.31.5") || (ver >= "1.30.5" && ver < "1.31.0");

  # If not, the package needs to be attached.
  if (!properS3methods) .requirePkg("aroma.light", quietly=TRUE);

  # Done?
  if (is.null(name)) return(invisible(NULL))


  # Retrieve a particular function?
  ns <- getNamespace("aroma.light");
  res <- get(name, mode=mode, envir=ns);

  # Function-specific patches?
  if (name == "normalizeTumorBoost") {
    # Workaround for older versions of aroma.light::normalizeTumorBoost()
    # assuming that the 'R.utils' package is attached. /HB 2013-09-20
    if (!properS3methods) .requirePkg("R.utils", quietly=TRUE);
  }

  invisible(res);
} # .useAromaLight()


.onAttach <- function(libname, pkgname) {
  pkg <- Package(pkgname);

  # Copy some pre-memoized CBS-parameter calculations to the 'R.cache'
  # cache.  This speeds up the calculation for common CBS use cases.
  .prememoize();

  # Inform user if DNAcopy is missing
  if (isPackageInstalled("DNAcopy")) {
    .requirePkg("DNAcopy");
  } else {
    msg <- "The Bioconductor package 'DNAcopy' is not installed. Please see http://www.bioconductor.org/ on how to install it, or try calling PSCBS::installDNAcopy().";
    hrule <- paste(rep("*", times=getOption("width", 80L)-1L), collapse="");
    packageStartupMessage(sprintf("%s\nNOTE: %s\n%s", hrule, msg, hrule));
  }

  startupMessage(pkg);
}


############################################################################
# HISTORY:
# 2013-09-27
# o Added .useAromaLight() to simplify backward compatibility.
# o Added .requirePkg() from the R.rsp package.
# 2011-07-23
# o Added a namespace to the package.
############################################################################
