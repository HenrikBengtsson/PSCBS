# aroma.core

<details>

* Version: 3.2.2
* GitHub: https://github.com/HenrikBengtsson/aroma.core
* Source code: https://github.com/cran/aroma.core
* Date/Publication: 2021-01-05 05:10:12 UTC
* Number of recursive dependencies: 47

Run `revdep_details(, "aroma.core")` for more info

</details>

## In both

*   checking package dependencies ... NOTE
    ```
    Packages suggested but not available for checking:
      'sfit', 'expectile', 'HaarSeg', 'mpcbs'
    ```

# PureCN

<details>

* Version: 1.20.0
* GitHub: https://github.com/lima1/PureCN
* Source code: https://github.com/cran/PureCN
* Date/Publication: 2020-10-27
* Number of recursive dependencies: 148

Run `revdep_details(, "PureCN")` for more info

</details>

## In both

*   checking tests ...
    ```
     ERROR
    Running the tests in ‘tests/testthat.R’ failed.
    Last 50 lines of output:
        1. └─PureCN::annotateTargets(...) test_annotateTargets.R:17:4
        2.   └─PureCN:::.checkSeqlevelStyle(x, txdb, "txdb", "interval file")
        3.     ├─GenomeInfoDb::`seqlevelsStyle<-`(`*tmp*`, value = refSeqlevelStyle[1])
        4.     └─GenomeInfoDb::`seqlevelsStyle<-`(`*tmp*`, value = refSeqlevelStyle[1])
        5.       ├─GenomeInfoDb::`seqlevelsStyle<-`(`*tmp*`, value = value)
        6.       └─GenomeInfoDb::`seqlevelsStyle<-`(`*tmp*`, value = value)
        7.         └─BiocGenerics::mapply(...)
    ...
       17.                   ├─base::do.call(...)
       18.                   └─(function (UCSC_chrom_info, assembly_accession, AssemblyUnits = NULL, ...
       19.                     └─GenomeInfoDb::getChromInfoFromNCBI(assembly_accession, assembly.units = AssemblyUnits)
       20.                       └─GenomeInfoDb:::.get_NCBI_chrom_info_from_accession(...)
       21.                         └─GenomeInfoDb:::fetch_assembly_report(accession)
       22.                           └─GenomeInfoDb:::.make_assembly_report_URL(assembly_accession)
      
      [ FAIL 2 | WARN 11 | SKIP 1 | PASS 340 ]
      Error: Test failures
      Execution halted
    ```

*   checking package dependencies ... NOTE
    ```
    Package which this enhances but not available for checking: ‘genomicsdb’
    ```

*   checking installed package size ... NOTE
    ```
      installed size is  8.7Mb
      sub-directories of 1Mb or more:
        doc       3.5Mb
        extdata   3.8Mb
    ```

