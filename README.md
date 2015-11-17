# PSCBS: Analysis of Parent-Specific DNA Copy Numbers


## Parallel processing
The package supports segmentation of the chromosomes in parallel
(asynchronously) via [futures](https://cran.r-project.org/package=future)
by adding the following
```r
future::plan("multicore")
```
to the beginning of the PSCBS script.  Everything else will work the
same.  To reset to non-parallel (synchronously) processing, use
`future::plan("eager")`.

An alternative to editing the R script is to set environment variable
`R_FUTURE_PLAN`, e.g.
```sh
export R_FUTURE_PLAN=multicore
```
To control the maximum number of cores the multicore processing may
use set environment variable `MC_CORES`, e.g.
```sh
export MC_CORES=4
```
This variable is defined and read by the 'parallel' package when it
is loaded (so not when R itself is started) and used to set options
`mc.cores`, which is acknowledged by `future::plan("multicore")`.

This above also be set in the cross-platform `~/.Renviron` file as:
```r
R_FUTURE_PLAN=multicore
MC_CORES=4
```



## Installation
R package PSCBS is available on [CRAN](http://cran.r-project.org/package=PSCBS) and can be installed in R as:
```r
install.packages('PSCBS')
```

### Pre-release version

To install the pre-release version that is available in branch `develop`, use:
```r
source('http://callr.org/install#HenrikBengtsson/PSCBS@develop')
```
This will install the package from source.  


## Software status

| Resource:     | CRAN        | Travis CI     | Appveyor         |
| ------------- | ------------------- | ------------- | ---------------- |
| _Platforms:_  | _Multiple_          | _Linux_       | _Windows_        |
| R CMD check   | <a href="http://cran.r-project.org/web/checks/check_results_PSCBS.html"><img border="0" src="http://www.r-pkg.org/badges/version/PSCBS" alt="CRAN version"></a> | <a href="https://travis-ci.org/HenrikBengtsson/PSCBS"><img src="https://travis-ci.org/HenrikBengtsson/PSCBS.svg" alt="Build status"></a> | <a href="https://ci.appveyor.com/project/HenrikBengtsson/pscbs"><img src="https://ci.appveyor.com/api/projects/status/github/HenrikBengtsson/PSCBS?svg=true" alt="Build status"></a> |
| Test coverage |                     | <a href="https://coveralls.io/r/HenrikBengtsson/PSCBS"><img src="https://coveralls.io/repos/HenrikBengtsson/PSCBS/badge.svg?branch=develop" alt="Coverage Status"/></a>   |                  |
