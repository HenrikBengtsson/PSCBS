library("PSCBS")
okind <- RNGkind()[1L]

sample1 <- function() { sample(0:9, size=1L) }

genv <- globalenv()

if (!is.null(genv$.Random.seed)) rm(list=".Random.seed", envir=genv, inherits=FALSE)

seed0 <- genv$.Random.seed
stopifnot(is.null(seed0))

## Get random seed
seed <- randomSeed("get")
stopifnot(identical(seed, seed0))

## Repeat after new sample
y1 <- sample1()
message(sprintf("Random number: %d", y1))
seed1 <- randomSeed("get")
stopifnot(!identical(seed1, seed0))

## Fixed random seed
randomSeed("set", seed=42L)
seed2 <- randomSeed("get")
stopifnot(!identical(seed2, seed1))

y2 <- sample1()
message(sprintf("Random number: %d (with random seed = 42L)", y2))

## Reset to previous state
randomSeed("reset")
seed3 <- randomSeed("get")
stopifnot(identical(seed3, seed1))

## Reset to NULL
randomSeed("set", seed=NULL)
seed4 <- randomSeed("get")
stopifnot(is.null(seed4))

y3 <- sample1()
message(sprintf("Random number: %d", y3))

## Fixed random seed
randomSeed("set", seed=42L)
y4 <- sample1()
message(sprintf("Random number: %d (with random seed = 42L)", y4))
stopifnot(identical(y4, y2))



## L'Ecuyer-CMRG: Random number generation for parallel processing
RNGkind("L'Ecuyer-CMRG")

randomSeed("set", seed=42L)
seed0 <- randomSeed("get")
seeds0 <- lapply(1:10, FUN=function(i) randomSeed("advance"))

## Assert reproducible .Random.seed stream
randomSeed("set", seed=42L)
seed1 <- randomSeed("get")
seeds1 <- lapply(1:10, FUN=function(i) randomSeed("advance"))
stopifnot(identical(seed1, seed0))
stopifnot(identical(seeds1, seeds0))


randomSeed("set", seed=42L)
y0 <- sapply(1:10, FUN=function(ii) {
  randomSeed("advance")
  sample1()
})
print(y0)
randomSeed("reset")

randomSeed("set", seed=42L)
y1 <- sapply(1:10, FUN=function(ii) {
  randomSeed("advance")
  sample1()
})
print(y1)
stopifnot(identical(y1, y0))
randomSeed("reset")


## Cleanup
RNGkind(okind)
