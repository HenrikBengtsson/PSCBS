library("PSCBS")
library("future")
oplan <- plan()

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

strategies <- c("eager", "lazy", "multicore")

## AGILE: Lazy futures fails with globals (<= 0.4.0)
if (packageVersion("globals") <= "0.4.0") {
  strategies <- setdiff(strategies, "lazy")
}

pkg <- "async"
if (R.utils::isPackageInstalled(pkg)) {
  library(pkg, character.only=TRUE)
  backend("local")
  strategies <- c(strategies, "batchjobs")
}

message("Future strategies: ", paste(sQuote(strategies), collapse=", "))

for (strategy in strategies) {
  message(sprintf("*** segmentByCBS() using '%s' futures ...", strategy))
  plan(strategy)
  fit <- segmentByCBS(data, verbose=TRUE)
  fits[[strategy]] <- fit
  stopifnot(all.equal(fit, fits[[1]]))
}

## Cleanup
plan(oplan)
