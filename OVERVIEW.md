## Parallel processing

The package supports segmentation of the chromosomes in parallel
(asynchronously) via [futures](https://cran.r-project.org/package=future)
by adding the following

```r
future::plan("multisession")
```

to the beginning of the PSCBS script.  Everything else will work the
same.  To reset to non-parallel processing, use `future::plan("sequential")`.

To configure this automatically whenever the package is loaded, see
future vignette '[A Future for R: Controlling Default Future Strategy](https://cran.r-project.org/web/packages/future/vignettes/future-5-startup.html)'.
