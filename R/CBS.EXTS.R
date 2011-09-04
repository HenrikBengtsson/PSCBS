###########################################################################/**
# @set class=DNAcopy
# @RdocMethod as.CBS
#
# @title "Coerces a DNAcopy object to a CBS object"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{fit}{A @see "DNAcopy" object (of the \pkg{DNAcopy} package.)}
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns a @see "CBS" object.
# }
#
# @author
#
# \seealso{
#   \code{\link[PSCBS:as.DNAcopy.CBS]{as.DNAcopy()}}.
#   @seeclass
# }
#
# @keyword internal
#*/########################################################################### 
setMethodS3("as.CBS", "DNAcopy", function(fit, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup the 'data' field
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  data <- fit$data;

  # Rename column names
  names <- colnames(data);
  names[names == "chrom"] <- "chromosome";
  names[names == "maploc"] <- "x";
  colnames(data) <- names;

  # Drop 'CNA' class and DNAcopy attributes
  class(data) <- c("data.frame");
  attr(data, "data.type") <- NULL;


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup the 'output' field
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  output <- fit$output;

  # Rename column names
  names <- colnames(output);
  names[names == "ID"] <- "id";
  names[names == "chrom"] <- "chromosome";
  names[names == "loc.start"] <- "start";
  names[names == "loc.end"] <- "end";
  names[names == "num.mark"] <- "nbrOfLoci";
  names[names == "seg.mean"] <- "mean";
  colnames(output) <- names;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup up 'CBS' object
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  res <- list();
  res$data <- data;
  res$output <- output;
  res$params <- list();
  class(res) <- c("CBS");

  res;
}) # as.CBS()



setMethodS3("extractByChromosomes", "CBS", function(x, chromosomes, ...) {
  # To please R CMD check
  this <- x;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'chromosomes':
  chromosomes <- Arguments$getIntegers(chromosomes, disallow=c("NaN", "Inf"));
  stopifnot(all(is.element(chromosomes, getChromosomes(this))));

  # Always extract in order
  chromosomes <- unique(chromosomes);
  chromosomes <- sort(chromosomes);

  # Allocate results
  res <- this;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Locus data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  chromosome <- NULL; rm(chromosome); # To please R CMD check
  data <- this$data;
  class <- class(data);
  class(data) <- "data.frame";
  data <- subset(data, chromosome %in% chromosomes);
  class(data) <- class;
  res$data <- data;


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Segmentation data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify rows to subset
  rows <- which(is.element(res$output$chromosome, chromosomes));
  for (field in c("output", "segRows")) {
    res[[field]] <- res[[field]][rows,,drop=FALSE];
  }

  # Identify chromosome offsets
  chrStarts <- match(getChromosomes(this), this$data$chromosome);
  chrEnds <- c(chrStarts[-1]-1L, nrow(this$data));
  chrLengths <- chrEnds - chrStarts + 1L;

  chrLengthsExcl <- chrLengths;

  keep <- match(chromosomes, getChromosomes(this));
  chrLengthsExcl[keep] <- 0L;
  cumChrLengthsExcl <- cumsum(chrLengthsExcl);

  shifts <- cumChrLengthsExcl[keep];
  stopifnot(all(is.finite(shifts)));

  # Adjust indices
  for (cc in seq(along=chromosomes)) {
    chromosome <- chromosomes[cc];
    shift <- shifts[cc];
    # Nothing to do?
    if (shift == 0) next;
    for (field in c("segRows")) {
      segRows <- res[[field]];
      rows <- which(res$output$chromosome == chromosome);
      segRows[rows,] <- segRows[rows,] - shift;
      res[[field]] <- segRows;
    }
  }

  res;
})


setMethodS3("subset", "CBS", function(x, chromlist=NULL, ...) {
  extractByChromosomes(x, chromosomes=chromlist, ...);
}, private=TRUE)


############################################################################
# HISTORY:
# 2011-09-03
# o Added as.CBS() for DNAcopy to coerce a DNAcopy object to a CBS object.
# 2011-09-02
# o Added extractByChromosomes() for CBS.
# o Added subset() for CBS for backward compatibility.
# o Added nbrOfLoci(), getChromosomes() and getSampleNames() for CBS.
# 2010-11-19
# o Added append() for CBS objects.
############################################################################
