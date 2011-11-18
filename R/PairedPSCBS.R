###########################################################################/**
# @RdocClass PairedPSCBS
#
# @title "The PairedPSCBS class"
#
# \description{
#  @classhierarchy
#
#  A PairedPSCBS is an object containing the results from the
#  Paired PSCBS method.
# }
# 
# \usage{PairedPSCBS(fit=list(), ...)}
#
# \arguments{
#   \item{fit}{A @list structure containing the Paired PSCBS results.}
#   \item{...}{Not used.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
# 
# @author
#
# \seealso{
#   The @see "segmentByPairedPSCBS" method returns an object of this class.
# }
#*/###########################################################################
setConstructorS3("PairedPSCBS", function(fit=list(), ...) {
  # Argument 'fit':
  if (!is.list(fit)) {
    throw("Argument 'fit' is not a list: ", class(fit)[1]);
  }

  extend(PSCBS(fit=fit, ...), "PairedPSCBS");
})


setMethodS3("getSegmentSizes", "PairedPSCBS", function(fit, by=c("length", "count"), ...) {
  by <- match.arg(by);

  data <- getSegments(fit, ...);
  if (by == "length") {
    res <- data[["tcnEnd"]]-data[["tcnStart"]]+1L;
  } else if (by == "count") {
    res <- data[["tcnNbrOfLoci"]];
  }
  res;
})


setMethodS3("updateMeans", "PairedPSCBS", function(fit, from=c("loci", "segments"), adjustFor=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'from':
  from <- match.arg(from);

  # Argument 'adjustFor':
  if (!is.null(adjustFor)) {
    adjustFor <- Arguments$getCharacters(adjustFor);
    adjustFor <- tolower(adjustFor);
    knownValues <- c("ab", "loh", "roh");
    adjustFor <- match.arg(adjustFor, choices=knownValues, several.ok=TRUE);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }
 
  verbose && enter(verbose, "Updating mean level estimates");
  verbose && cat(verbose, "Adjusting for:");
  verbose && print(verbose, adjustFor);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Extract the segmentation results
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  segs <- getSegments(fit, splitters=TRUE);
  nbrOfSegments <- nrow(segs);
  verbose && cat(verbose, "Number of segments: ", nbrOfSegments);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Assert that adjustments can be made
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  if (is.element("ab", adjustFor)) {
    if (!is.element("abCall", names(segs))) {
      adjustFor <- setdiff(adjustFor, "ab");
      throw("Cannot adjust for AB, because they haven't been called.");
    }
  }

  if (is.element("loh", adjustFor)) {
    if (!is.element("lohCall", names(segs))) {
      adjustFor <- setdiff(adjustFor, "loh");
      throw("Cannot adjust for LOH, because they haven't been called.");
    }
  }

  if (is.element("roh", adjustFor)) {
    if (!is.element("rohCall", names(segs))) {
      adjustFor <- setdiff(adjustFor, "roh");
      throw("Cannot adjust for ROH, because they haven't been called.");
    }
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Update the (TCN,DH) mean levels from locus-level data?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  if (from == "loci") {
    data <- getLocusData(fit);
    chromosome <- data$chromosome;
    x <- data$x;
    CT <- data$CT;
    rho <- data$rho;
    muN <- data$muN;
    isSnp <- !is.na(muN);
    isHet <- isSnp & (muN == 1/2);

    isSplitter <- isSegmentSplitter(fit);
    for (ss in seq(length=nbrOfSegments)[!isSplitter]) {
      verbose && enter(verbose, sprintf("Segment %d of %d", ss, nbrOfSegments));
      seg <- segs[ss,];
      verbose && print(verbose, seg);
  
      chr <- seg[["chromosome"]];
      chrTag <- sprintf("chr%02d", chr);
  
      for (what in c("tcn", "dh")) {
        keys <- paste(what, c("Start", "End"), sep="");
        xStart <- seg[[keys[1]]];
        xEnd <- seg[[keys[2]]];
        verbose && printf(verbose, "[xStart,xEnd] = [%.0f,%.0f]\n", xStart, xEnd);
        # Nothing todo?
        if (is.na(xStart) && is.na(xEnd)) {
          next;
        }
  
        stopifnot(xStart <= xEnd);
  
        # (b) Identify units
        units <- which(chromosome == chr & xStart <= x & x <= xEnd);
  
        # (c) Adjust for missing values
        if (what == "tcn") {
          value <- CT;
        } else if (what == "dh") {
          value <- rho;
        }
        keep <- which(!is.na(value[units]));
        units <- units[keep];
  
        # (d) Update mean
        gamma <- mean(value[units]);
  
        # Sanity check
        stopifnot(length(units) == 0 || !is.na(gamma));
  
        # Update the segment boundaries, estimates and counts
        key <- paste(what, "Mean", sep="");
        seg[[key]] <- gamma;
      }
  
      verbose && print(verbose, seg);
  
      segs[ss,] <- seg;
  
      verbose && exit(verbose);
    } # for (ss ...)
  } # if (from ...)



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Adjust segment means from various types of calls
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  if (length(adjustFor) > 0) {
    verbose && enter(verbose, "Adjusting segment means");
    verbose && cat(verbose, "Adjusting for:");
    verbose && print(verbose, adjustFor);

    if (is.element("ab", adjustFor)) {
      verbose && enter(verbose, "Adjusting for AB");
      calls <- segs$abCall;
      segs$dhMean[calls] <- 1/2;
      verbose && exit(verbose);
    }

    if (is.element("loh", adjustFor)) {
      verbose && enter(verbose, "Adjusting for LOH");
      calls <- segs$lohCall;
      segs$dhMean[calls] <- 0;
      verbose && exit(verbose);
    }

    if (is.element("roh", adjustFor)) {
      verbose && enter(verbose, "Adjusting for ROH");
      calls <- segs$rohCall;
      segs$dhMean[calls] <- as.double(NA);
      verbose && exit(verbose);
    }

    verbose && exit(verbose);
  } # if (length(adjustFor) > 0)



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Update (C1,C2) mean levels
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  verbose && enter(verbose, "Update (C1,C2) per segment");
  # Append (C1,C2) estimates
  tcn <- segs$tcnMean;
  dh <- segs$dhMean;
  C1 <- 1/2*(1-dh)*tcn;
  C2 <- tcn - C1;
  segs$c1Mean <- C1;
  segs$c2Mean <- C2;
  verbose && exit(verbose);


  # Return results
  res <- fit;
  res$output <- segs;

  verbose && exit(verbose);

  res;
}, private=TRUE) # updateMeans()




setMethodS3("resegment", "PairedPSCBS", function(fit, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Resegmenting a ", class(fit)[1], " object");

  # Use the locus-level data of the PairedPSCBS object
  data <- getLocusData(fit);
  class(data) <- "data.frame";
  drop <- c("rho", "betaTN", "index");
  keep <- !is.element(colnames(data), drop);
  data <- data[,keep];
  verbose && str(verbose, data);

  verbose && cat(verbose, "Number of loci: ", nrow(data));


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup arguments to be passed
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Overriding default arguments");
  segFcnName <- "segmentByCBS";
  segFcn <- getMethodS3(segFcnName, "default");

  # (a) The default arguments
  formals <- formals(segFcn);

  formals <- formals[!sapply(formals, FUN=is.language)];
  formals <- formals[!sapply(formals, FUN=is.name)];
  drop <- c("chromosome", "x", "w", "CT", "betaT", "betaN", "muN", "...");
  keep <- !is.element(names(formals), drop);
  formals <- formals[keep];

  # (b) The arguments used in previous fit
  params <- fit$params;
  keep <- is.element(names(params), names(formals));
  params <- params[keep];
  # Don't trust 'tbn'!  TODO. /HB 20111117
  params$tbn <- NULL;

  # (c) The arguments in '...'
  userArgs <- list(..., verbose=verbose);

  # (d) Merge
  args <- formals;  
  args2 <- append(params, userArgs);
  for (kk in seq(along=args2)) {
    value <- args2[[kk]];
    if (!is.null(value)) {
      key <- names(args2)[kk];
      if (!is.null(key)) {
        args[[key]] <- value;
      } else {
        args <- append(args, list(value));
      }
    }
  } # for (key ...)
  verbose && str(verbose, args[names(args) != "verbose"]);

  args <- append(list(data), args);
  verbose && cat(verbose, "Arguments with data:");
  verbose && str(verbose, args[names(args) != "verbose"]);
  verbose && exit(verbose);

  verbose && enter(verbose, sprintf("Calling %s()", segFcnName));
  fit <- do.call(segFcnName, args);
  verbose && exit(verbose);

  verbose && exit(verbose);

  fit;
}, protected=TRUE) # resegment()



##############################################################################
# HISTORY
# 2011-11-17
# o Added resegment() for PairedPSCBS for easy resegmentation.
# 2011-11-12
# o Added arguments 'from' and 'adjustFor' to updateMeans().
# 2011-01-16
# o BUG FIX: updateMeans() save to the incorrect column names.
# 2011-01-12
# o Added updateMeans() for PairedPSCBS.
# 2011-10-02
# o CLEANUP: Moved print() and as.data.frame() to PSCBS.
# o Added Rdoc help.
# o Now the constructor of PairedPSCBS calls that of PSCBS.
# 2011-06-28
# o DOCUMENTATION: Added Rd help for as.data.frame() of PairedPSCBS.
# 2011-04-08
# o Added formal constructor for the PairedPSCBS class.
# o Created.
##############################################################################
