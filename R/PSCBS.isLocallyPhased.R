setMethodS3("isLocallyPhased", "PSCBS", function(fit, ...) {
  segs <- getSegments(fit)
  is.element("c1c2Swap", names(segs))
})


##############################################################################
# HISTORY
# 2014-03-24 [HB]
# o Added isLocallyPhased() for PSCBS.
# o Created.
##############################################################################
