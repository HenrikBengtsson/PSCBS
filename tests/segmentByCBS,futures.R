library("PSCBS")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Simulating copy-number data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
set.seed(0xBEEF)

# Number of loci
J <- 1000

mu <- double(J)
mu[200:300] <- mu[200:300] + 1
mu[350:400] <- NA # centromere
mu[650:800] <- mu[650:800] - 1
eps <- rnorm(J, sd=1/2)
y <- mu + eps
x <- sort(runif(length(y), max=length(y))) * 1e5
w <- runif(J)
w[650:800] <- 0.001

## Create multiple chromosomes
data <- list()
for (cc in 1:3) {
  data[[cc]] <- data.frame(chromosome=cc, y=y, x=x)
}
data <- Reduce(rbind, data)
str(data)

## Segment
fits <- list()


message("*** segmentByCBS() via futures ...")

message("*** segmentByCBS() via 'lazy' futures without attaching 'future'...")
future::plan("lazy")
print(future::plan)
fitL <- segmentByCBS(data, verbose=TRUE)
print(fitL)


message("*** segmentByCBS() via futures with 'future' attached...")
library("future")
oplan <- plan()

strategies <- c("eager", "lazy", "multicore")

## Also test BatchJobs futures in async?
pkg <- "async"
if (R.utils::isPackageInstalled(pkg)) {
  library(pkg, character.only=TRUE)
  backend("local")
  strategies <- c(strategies, "batchjobs")
}

message("Future strategies to test: ", paste(sQuote(strategies), collapse=", "))

for (strategy in strategies) {
  message(sprintf("- segmentByCBS() using '%s' futures ...", strategy))
  plan(strategy)
  fit <- segmentByCBS(data, verbose=TRUE)
  fits[[strategy]] <- fit
  stopifnot(all.equal(fit, fits[[1]]))
  stopifnot(all.equal(fit, fitL))
}

## Cleanup
plan(oplan)
