# CRAN submission PSCBS 0.45.0
on 2015-09-11

Changes related to R/CRAN updates:

* Explicitly importing core R functions.

Thank you in advance


## Notes not sent to CRAN
PSCBS have been verified using `R CMD check --as-cran` on:

* Platform x86_64-pc-linux-gnu (64-bit):
  - R version 3.1.1 (2014-07-10)
  - R version 3.1.3 (2015-03-09)
  - R version 3.2.2 (2015-08-14)
  - R version 3.2.2 Patched (2015-09-09 r69342)
* Platform x86_64-w64-mingw32/x64 (64-bit):
  - R version 3.1.3 (2015-03-09)
  - R version 3.2.2 Patched (2015-09-08 r69332)
  - R Under development (unstable) (2015-09-07 r69315)

It has also verified using the <http://win-builder.r-project.org/> service.

Moreover, the updates cause no issues for any of the following 2
reverse dependencies on CRAN and Bioconductor, which have been tested
with `R CMD check --as-cran`: aroma.cn 1.6.0 and aroma.core 2.13.1.
