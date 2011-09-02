setMethodS3("nbrOfSegments", "DNAcopy", function(fit, ...) {
  segs <- fit$output;
  nrow(segs);
})

setMethodS3("nbrOfLoci", "DNAcopy", function(fit, ...) {
  nrow(fit$data);
})

setMethodS3("nbrOfSamples", "DNAcopy", function(fit, ...) {
  length(getSampleNames(fit, ...));
})

setMethodS3("getSampleNames", "DNAcopy", function(fit, ...) {
  names <- colnames(fit$data);
  names <- setdiff(names, c("chrom", "maploc"));
  names;
})

setMethodS3("getChromosomes", "DNAcopy", function(fit, ...) {
  chromosomes <- fit$data$chrom;
  sort(unique(chromosomes));
})


############################################################################
# HISTORY:
# 2011-09-02
# o ROBUSTNESS: Now getSampleNames() drops columns 'chrom' and 'maploc',
#   instead of assuming their positions.
# o ROBUSTNESS: Now nbrOfSamples() utilizes getSampleNames().
# o Added nbrOfSegments(), nbrOfLoci(), nbrOfSamples(), getSampleNames()
#   and getChromosomes() for DNAcopy.
# HISTORY FROM PRIVATE SCRIPTS:
# 2011-07-18
# o Added getSampleNames().
# o Added plotTracks() for DNAcopy.
# o Added nbrOfSegments(), nbrOfLoci() and nbrOfSamples().
############################################################################
