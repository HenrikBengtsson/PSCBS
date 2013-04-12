verbose <- Arguments$getVerbose(-10*interactive(), timestamp=TRUE)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Load SNP microarray data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
data <- PSCBS::exampleData("paired.chr01")
str(data)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Paired PSCBS segmentation
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Drop single-locus outliers
dataS <- dropSegmentationOutliers(data)

# Speed up example by segmenting fewer loci
#dataS <- dataS[seq(from=1, to=nrow(data), by=5),]

str(dataS)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Paired PSCBS segmentation
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
gaps <- findLargeGaps(dataS, minLength=2e6)
knownSegments <- gapsToSegments(gaps, dropGaps=TRUE)

fit <- segmentByPairedPSCBS(dataS, knownSegments=knownSegments,
                            seed=0xBEEF, verbose=verbose)
print(fit)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Paired PSCN mean levels of hets and homs separately
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Fix segments as above
segs <- getSegments(fit)[,c("chromosome", "tcnStart", "tcnEnd")]
names(segs)[2:3] <- c("start", "end")

# Identify heterozygous SNPs
isHet <- (getLocusData(fit)$muN == 1/2)
 
fitList <- list()
for (key in c("homozygotes", "heterozygotes")) {
  dataT <- dataS
  excl <- if(key == "homozygotes") isHet else !isHet
  dataT$CT[excl] <- NA

  fitT <- segmentByPairedPSCBS(dataT, knownSegments=segs,
                               undoTCN=Inf, undoDH=Inf, verbose=verbose)
  fitList[[key]] <- fitT
} # for (key ...)


# Plot results
subplots(2, ncol=1)
for (key in names(fitList)) {
  fitT <- fitList[[key]]
  plotTracks(fitT, tracks="tcn,c1,c2", add=FALSE)
  abline(h=0:6, lty=3, col="gray")
  title(main=key)
}
