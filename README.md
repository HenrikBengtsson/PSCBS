# PSCBS: Analysis of Parent-Specific DNA Copy Numbers


## Parallel processing
The package supports segmentation of the chromosomes in parallel
(asynchronously) via [futures](https://cran.r-project.org/package=future)
by adding the following
```r
future::plan("multiprocess")
```
to the beginning of the PSCBS script.  Everything else will work the
same.  To reset to non-parallel (synchronously) processing, use
`future::plan("eager")`.

To configure this automatically whenever the package is loaded, see
future vignette '[A Future for R: Controlling Default Future Strategy](https://cran.r-project.org/web/packages/future/vignettes/future-4-startup.html)'.


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

| Resource:     | CRAN        | Travis CI       | Appveyor         |
| ------------- | ------------------- | --------------- | ---------------- |
| _Platforms:_  | _Multiple_          | _Linux & macOS_ | _Windows_        |
| R CMD check   | <a href="https://cran.r-project.org/web/checks/check_results_PSCBS.html"><img border="0" src="http://www.r-pkg.org/badges/version/PSCBS" alt="CRAN version"></a> | <a href="https://travis-ci.org/HenrikBengtsson/PSCBS"><img src="https://travis-ci.org/HenrikBengtsson/PSCBS.svg" alt="Build status"></a>   | <a href="https://ci.appveyor.com/project/HenrikBengtsson/pscbs"><img src="https://ci.appveyor.com/api/projects/status/github/HenrikBengtsson/PSCBS?svg=true" alt="Build status"></a> |
| Test coverage |                     | <a href="https://codecov.io/gh/HenrikBengtsson/PSCBS"><img src="https://codecov.io/gh/HenrikBengtsson/PSCBS/branch/develop/graph/badge.svg" alt="Coverage Status"/></a>     |                  |
