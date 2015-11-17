# eCRAN submission PSCBS 0.60.0
on 2015-11-17

Version jump in minor is intentional.

Thanks in advance


## Notes not sent to CRAN
The package has been verified using `R CMD check --as-cran` on:

* Platform x86_64-pc-linux-gnu (64-bit):
  - R version 3.1.1 (2014-07-10)
  - R version 3.1.3 (2015-03-09)
  - R version 3.2.2 (2015-08-14)
  - R version 3.2.2 Patched (2015-11-13 r69635)
  - R Under development (unstable) (2015-10-13 r69635)

* Platform: x86_64-apple-darwin13.4.0 (64-bit):
  - R version 3.2.2 Patched (2015-09-05 r69301)
  
* Platform x86_64-w64-mingw32/x64 (64-bit):
  - R version 3.1.3 (2015-03-09)
  - R version 3.2.2 Patched (2015-11-13 r69636)
  - R Under development (unstable) (2015-11-16 r69640)

It has also verified using the <http://win-builder.r-project.org/> service.

Moreover, the updates cause no issues for any of the following 2
reverse dependencies on CRAN and Bioconductor, which have been tested
with `R CMD check --as-cran`: aroma.cn 1.6.1 and aroma.core 2.14.0.
