library("PSCBS")
rngSeeder <- PSCBS:::rngSeeder

sample1 <- function() { sample(0:9, size=1L) }

genv <- globalenv()

if (!is.null(genv$.Random.seed)) rm(list=".Random.seed", envir=genv, inherits=FALSE)

seed0 <- genv$.Random.seed

## Get random seed
seed <- rngSeeder("get")
stopifnot(identical(seed, seed0))

## Repeat after new sample
y1 <- sample1()
message(sprintf("Random number: %d", y1))
seed1 <- rngSeeder("get")
stopifnot(!identical(seed1, seed0))

## Fixed random seed
rngSeeder("set", seed=42L)
seed2 <- rngSeeder("get")
stopifnot(!identical(seed2, seed1))

y2 <- sample1()
message(sprintf("Random number: %d (with random seed = 42L)", y2))

## Reset to previous state
rngSeeder("reset")
seed3 <- rngSeeder("get")
stopifnot(!identical(seed3, seed1))

## Reset to NULL
rngSeeder("set", seed=NULL)
seed4 <- rngSeeder("get")
stopifnot(identical(seed4, seed0))

y3 <- sample1()
message(sprintf("Random number: %d", y3))

## Fixed random seed
rngSeeder("set", seed=42L)
y4 <- sample1()
message(sprintf("Random number: %d (with random seed = 42L)", y4))
stopifnot(identical(y4, y2))
