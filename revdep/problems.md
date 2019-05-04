# aroma.cn

<details>

* Version: 1.6.1
* Source code: https://github.com/cran/aroma.cn
* URL: http://www.aroma-project.org/, https://github.com/HenrikBengtsson/aroma.cn
* BugReports: https://github.com/HenrikBengtsson/aroma.cn/issues
* Date/Publication: 2015-10-28 00:08:16
* Number of recursive dependencies: 23

Run `revdep_details(,"aroma.cn")` for more info

</details>

## In both

*   checking whether package ‘aroma.cn’ can be installed ... NOTE
    ```
    Found the following notes/warnings:
      Non-staged installation was used
    See ‘/home/hb/repositories/PSCBS/revdep/checks/aroma.cn/new/aroma.cn.Rcheck/00install.out’ for details.
    ```

# aroma.core

<details>

* Version: 3.1.3
* Source code: https://github.com/cran/aroma.core
* URL: https://github.com/HenrikBengtsson/aroma.core, http://www.aroma-project.org/
* BugReports: https://github.com/HenrikBengtsson/aroma.core/issues
* Date/Publication: 2018-05-03 13:41:54 UTC
* Number of recursive dependencies: 45

Run `revdep_details(,"aroma.core")` for more info

</details>

## In both

*   checking package dependencies ... NOTE
    ```
    Packages suggested but not available for checking:
      'sfit', 'expectile', 'HaarSeg', 'mpcbs'
    ```

*   checking whether package ‘aroma.core’ can be installed ... NOTE
    ```
    Found the following notes/warnings:
      Non-staged installation was used
    See ‘/home/hb/repositories/PSCBS/revdep/checks/aroma.core/new/aroma.core.Rcheck/00install.out’ for details.
    ```

# EstMix

<details>

* Version: 1.0.1
* Source code: https://github.com/cran/EstMix
* Date/Publication: 2018-09-13 04:20:02 UTC
* Number of recursive dependencies: 15

Run `revdep_details(,"EstMix")` for more info

</details>

## In both

*   checking examples ... ERROR
    ```
    ...
    > ## short example
    > ##
    > #########################################################
    > ## first load the data
    > BAF <- example_data$BAF
    > LRR <- example_data$LRR ## In practice, the orignal LRR should be devided by 0.55
    > chr <- example_data$chr
    > loc <- example_data$x
    > GT <- example_data$GT
    > gt = (GT=='BB')*2+(GT=='AB')*1.5+(GT=='AA')-1;gt[gt==(-1)]=NA
    > 
    > ## then perform segmentation
    > gaps = PSCBS::findLargeGaps(x=loc,minLength=5e6,chromosome=chr)
    > if(!is.null(gaps)) knownSegments = PSCBS::gapsToSegments(gaps)
    > p <- 0.0001
    > fit <- PSCBS::segmentByPairedPSCBS(CT=2*2^LRR,betaT=BAF,muN=gt,chrom=chr,
    + knownSegments=knownSegments,tbn=FALSE,x=loc,seed=1, alphaTCN=p*.9,alphaDH=p*.1)
    Error in is.na(tcnSegRowsKK[, 1]) || is.na(dhSegRowsKK[, 1]) : 
      'length(x) = 3 > 1' in coercion to 'logical(1)'
    Calls: <Anonymous> -> segmentByPairedPSCBS.default -> .stop_if_not
    Execution halted
    ```

*   checking whether package ‘EstMix’ can be installed ... NOTE
    ```
    Found the following notes/warnings:
      Non-staged installation was used
    See ‘/home/hb/repositories/PSCBS/revdep/checks/EstMix/new/EstMix.Rcheck/00install.out’ for details.
    ```

# jointseg

<details>

* Version: 1.0.2
* Source code: https://github.com/cran/jointseg
* URL: https://github.com/mpierrejean/jointseg
* BugReports: https://github.com/mpierrejean/jointseg/issues
* Date/Publication: 2019-01-11 12:30:03 UTC
* Number of recursive dependencies: 42

Run `revdep_details(,"jointseg")` for more info

</details>

## In both

*   checking whether package ‘jointseg’ can be installed ... NOTE
    ```
    Found the following notes/warnings:
      Non-staged installation was used
    See ‘/home/hb/repositories/PSCBS/revdep/checks/jointseg/new/jointseg.Rcheck/00install.out’ for details.
    ```

# PureCN

<details>

* Version: 1.14.0
* Source code: https://github.com/cran/PureCN
* URL: https://github.com/lima1/PureCN
* Date/Publication: 2019-05-02
* Number of recursive dependencies: 120

Run `revdep_details(,"PureCN")` for more info

</details>

## In both

*   checking whether package ‘PureCN’ can be installed ... NOTE
    ```
    Found the following notes/warnings:
      Non-staged installation was used
    See ‘/home/hb/repositories/PSCBS/revdep/checks/PureCN/new/PureCN.Rcheck/00install.out’ for details.
    ```

*   checking installed package size ... NOTE
    ```
      installed size is  8.2Mb
      sub-directories of 1Mb or more:
        data      1.1Mb
        doc       2.7Mb
        extdata   3.8Mb
    ```

