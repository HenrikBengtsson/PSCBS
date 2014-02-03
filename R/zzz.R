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


# This function exists solely for the purpose of backward compatibility
# /HB 2014-02-04
.useDNAcopy <- function(name=NULL, mode="function") {
  if (!isPackageInstalled("DNAcopy")) {
    throw("Bioconductor package not installed: DNAcopy");
  }

  # Done?
  if (is.null(name)) return(invisible(NULL))

  # Retrieve a particular function?
  ns <- getNamespace("DNAcopy");
  res <- get(name, mode=mode, envir=ns);

  invisible(res);
} # .useDNAcopy()


.onLoad <- function(libname, pkgname) {
  ns <- getNamespace(pkgname);
  pkg <- Package(pkgname);
  # Assign '.PSCBS' object [since 'PSCBS' is a constructor/Class].
  name <- sprintf(".%s", pkgname);
  assign(name, pkg, envir=ns, inherits=FALSE);
} # .onLoad()


.onAttach <- function(libname, pkgname) {
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

  # Get '.PSCBS' object [since 'PSCBS' is a constructor/Class].
  name <- sprintf(".%s", pkgname);
  pkg <- get(name, envir=getNamespace(pkgname), inherits=FALSE);
  startupMessage(pkg);
}


############################################################################
# HISTORY:
# 2014-02-04
# o Added .useDNAcopy() to simplify backward compatibility.
# 2013-10-13
# o Added .onLoad() that creates Package '.PSCBS' object, which is
#   used in .onAttach().  This is a workaround for not allocating a
#   local Package on in .onAttach(), which then will be garbage
#   collected and finalize():d, which in turn can generate cyclic
#   loading of namespaces in R.oo (< 1.16.0).
# 2013-09-27
# o Added .useAromaLight() to simplify backward compatibility.
# o Added .requirePkg() from the R.rsp package.
# 2011-07-23
# o Added a namespace to the package.
############################################################################
