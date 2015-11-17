library("PSCBS")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Load SNP microarray data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
data <- PSCBS::exampleData("paired.chr01")


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Paired PSCBS segmentation
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Drop single-locus outliers
dataS <- dropSegmentationOutliers(data)

# Run light-weight tests by default
if (Sys.getenv("_R_CHECK_FULL_") == "") {
  # Use only every 5th data point
  dataS <- dataS[seq(from=1, to=nrow(data), by=5),]
  # Number of segments (for assertion)
  nSegs <- 4L
} else {
  # Full tests
  nSegs <- 11L
}

str(dataS)


## Create multiple chromosomes
data <- knownSegments <- list()
for (cc in 1:3) {
  dataS$chromosome <- cc
  data[[cc]] <- dataS
  knownSegments[[cc]] <- data.frame(
    chromosome = c(       cc, cc,        cc),
    start      = c(     -Inf, NA, 141510003),
    end        = c(120992603, NA,      +Inf)
  )
}
data <- Reduce(rbind, data)
str(data)
knownSegments <- Reduce(rbind, knownSegments)
str(knownSegments)


message("*** segmentByPairedPSCBS() via futures ...")

library("future")
oplan <- plan()

strategies <- c("eager", "lazy", "multicore")

## Test 'async' futures?
pkg <- "async"
if (require(pkg, character.only=TRUE)) {
  backend("local")
  strategies <- c(strategies, "batchjobs")
}

message("Future strategies to test: ", paste(sQuote(strategies), collapse=", "))

fits <- list()
for (strategy in strategies) {
  message(sprintf("- segmentByPairedPSCBS() using '%s' futures ...", strategy))
  plan(strategy)
  fit <- segmentByCBS(data, seed=0xBEEF, verbose=TRUE)
  fits[[strategy]] <- fit
  stopifnot(all.equal(fit, fits[[1]]))
}

message("*** segmentByPairedPSCBS() via futures ... DONE")

## Cleanup
plan(oplan)
rm(list=c("fits", "data", "fit"))
