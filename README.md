# PSCBS: Analysis of Parent-Specific DNA Copy Numbers


## Parallel processing
The package supports segmentation of the chromosomes in parallel
(asynchronously) via [futures](https://cran.r-project.org/package=future)
by adding the following
```r
future::plan("multiprocess")
```
to the beginning of the PSCBS script.  Everything else will work the
same.  To reset to non-parallel processing, use `future::plan("sequential")`.

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


## Contributions

This Git repository uses the [Git Flow](http://nvie.com/posts/a-successful-git-branching-model/) branching model (the [`git flow`](https://github.com/petervanderdoes/gitflow-avh) extension is useful for this).  The [`develop`](https://github.com/HenrikBengtsson/PSCBS/tree/develop) branch contains the latest contributions and other code that will appear in the next release, and the [`master`](https://github.com/HenrikBengtsson/PSCBS) branch contains the code of the latest release, which is exactly what is currently on [CRAN](https://cran.r-project.org/package=PSCBS).

Contributing to this package is easy.  Just send a [pull request](https://help.github.com/articles/using-pull-requests/).  When you send your PR, make sure `develop` is the destination branch on the [PSCBS repository](https://github.com/HenrikBengtsson/PSCBS).  Your PR should pass `R CMD check --as-cran`, which will also be checked by <a href="https://travis-ci.org/HenrikBengtsson/PSCBS">Travis CI</a> and <a href="https://ci.appveyor.com/project/HenrikBengtsson/pscbs">AppVeyor CI</a> when the PR is submitted.


## Software status

| Resource:     | CRAN        | Travis CI       | Appveyor         |
| ------------- | ------------------- | --------------- | ---------------- |
| _Platforms:_  | _Multiple_          | _Linux & macOS_ | _Windows_        |
| R CMD check   | <a href="https://cran.r-project.org/web/checks/check_results_PSCBS.html"><img border="0" src="http://www.r-pkg.org/badges/version/PSCBS" alt="CRAN version"></a> | <a href="https://travis-ci.org/HenrikBengtsson/PSCBS"><img src="https://travis-ci.org/HenrikBengtsson/PSCBS.svg" alt="Build status"></a>   | <a href="https://ci.appveyor.com/project/HenrikBengtsson/pscbs"><img src="https://ci.appveyor.com/api/projects/status/github/HenrikBengtsson/PSCBS?svg=true" alt="Build status"></a> |
| Test coverage |                     | <a href="https://codecov.io/gh/HenrikBengtsson/PSCBS"><img src="https://codecov.io/gh/HenrikBengtsson/PSCBS/branch/develop/graph/badge.svg" alt="Coverage Status"/></a>     |                  |
