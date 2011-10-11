setMethodS3("extractLocusLevelC1C2", "PairedPSCBS", function(fit, ...) {
  data <- getLocusData(fit);
  C <- data$CT;
  rho <- data$rho;
  C1 <- 1/2*(1-rho)*C;
  C2 <- C-C1;
  data.frame(C1=C1, C2=C2);
}, private=TRUE) # extractLocusLevelC1C2()


setMethodS3("extractLocusLevelTCN", "PairedPSCBS", function(fit, ...) {
  data <- getLocusData(fit);
  C <- data$CT;
}, private=TRUE) # extractLocusLevelTCN()



############################################################################
# HISTORY:
# 2011-10-8
# o ROBUSTIFICATION: Uses drop=FALSE in mergeTwoSegments() for PairedPSCBS.
# 2011-10-02
# o DOCUMENTATION: Added Rdoc help to mergeTwoSegments() & dropByRegions().
# o Added verbose statements to the above to functions.
# 2011-06-14
# o Updated code to recognize new column names.
# 2011-01-18
# o BUG FIX: Fields 'tcnSegRows' and 'dhSegRows' were not updated by
#   mergeTwoSegments() for PairedPSCBS.
# 2011-01-14
# o Moved extractByRegions() and estimateStdDevForHeterozygousBAF() to
#   psCBS v0.9.36.
# o Now extractByRegions() utilizes the 'segRows' field.
# o Added estimateStdDevForHeterozygousBAF().
# 2011-01-12
# o Added updateMeans() for PairedPSCBS.
# o Added dropByRegions().
# o Added extractByRegions() and extractByRegion().
# o Created.
############################################################################
