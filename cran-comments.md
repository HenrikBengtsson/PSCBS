# CRAN submission PSCBS 0.66.0

on 2021-10-22

Thanks in advance


## Notes not sent to CRAN

### R CMD check validation

The package has been verified using `R CMD check --as-cran` on:

| R version     | GitHub | R-hub      | mac/win-builder |
| ------------- | ------ | ---------- | --------------- |
| 3.6.x         | L      |            |                 |
| 4.0.x         | L      | L          |                 |
| 4.1.x         | L M W  | L M M1 S W | M1 W            |
| devel         | L M    | L        W |    W            |

*Legend: OS: L = Linux, S = Solaris, M = macOS, M1 = macOS M1, W = Windows*


R-hub checks:

```
> res <- rhub::check(platform = c(
  "debian-clang-devel", "debian-gcc-patched", "linux-x86_64-centos-epel",
  "macos-highsierra-release-cran", "macos-m1-bigsur-release",
  "solaris-x86-patched-ods", "windows-x86_64-devel", "windows-x86_64-release"))
> res

── PSCBS 0.66.0: OK

  Build ID:   PSCBS_0.66.0.tar.gz-fc019421033248ddad6126c3590edc2a
  Platform:   Debian Linux, R-devel, clang, ISO-8859-15 locale
  Submitted:  27m 57.1s ago
  Build time: 17m 28.2s

0 errors ✔ | 0 warnings ✔ | 0 notes ✔

── PSCBS 0.66.0: OK

  Build ID:   PSCBS_0.66.0.tar.gz-00739df3ce16490db4009bfa24bde017
  Platform:   Debian Linux, R-patched, GCC
  Submitted:  27m 57.1s ago
  Build time: 14m 36.2s

0 errors ✔ | 0 warnings ✔ | 0 notes ✔

── PSCBS 0.66.0: ERROR

  Build ID:   PSCBS_0.66.0.tar.gz-90ad4020be2a47f1996c7e649bbe2ce3
  Platform:   CentOS 8, stock R from EPEL
  Submitted:  27m 57.1s ago
  Build time: 10m 44.2s

❯ checking examples ... ERROR
  Running examples in ‘PSCBS-Ex.R’ failed
  The error most likely occurred in:
  
  > ### Name: segmentByNonPairedPSCBS
  > ### Title: Segment total copy numbers and allele B fractions using the
  > ###   Non-paired PSCBS method
  > ### Aliases: segmentByNonPairedPSCBS.default segmentByNonPairedPSCBS
  > ###   segmentByNonPairedPSCBS.data.frame
  > ###   segmentByNonPairedPSCBS.PairedPSCBS segmentByNonPairedPSCBS
  > ### Keywords: IO
  > 
  > ### ** Examples
  > 
  > verbose <- R.utils::Arguments$getVerbose(-10*interactive(), timestamp=TRUE)
  > 
  > # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  > # Load SNP microarray data
  > # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  > data <- PSCBS::exampleData("paired.chr01")
  > str(data)
  'data.frame':	73346 obs. of  6 variables:
   $ chromosome: int  1 1 1 1 1 1 1 1 1 1 ...
   $ x         : int  1145994 2224111 2319424 2543484 2926730 2941694 3084986 3155127 3292731 3695086 ...
   $ CT        : num  1.625 1.071 1.406 1.18 0.856 ...
   $ betaT     : num  0.757 0.771 0.834 0.778 0.229 ...
   $ CN        : num  2.36 2.13 2.59 1.93 1.71 ...
   $ betaN     : num  0.827 0.875 0.887 0.884 0.103 ...
  > 
  > 
  > # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  > # Paired PSCBS segmentation
  > # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  > # Drop single-locus outliers
  > dataS <- dropSegmentationOutliers(data)
  > 
  > # Speed up example by segmenting fewer loci
  > dataS <- dataS[seq(from=1, to=nrow(data), by=20),]
  > 
  > str(dataS)
  'data.frame':	3668 obs. of  6 variables:
   $ chromosome: int  1 1 1 1 1 1 1 1 1 1 ...
   $ x         : int  1145994 4276892 5034491 6266412 8418532 11211748 13928296 14370144 15014887 16589707 ...
   $ CT        : num  1.63 1.16 1.35 1.39 1.55 ...
   $ betaT     : num  0.7574 0.0576 0.8391 0.7917 0.8141 ...
   $ CN        : num  2.36 2.32 2.33 2.97 2.31 ...
   $ betaN     : num  0.8274 0.0421 0.9406 0.869 0.6045 ...
  > 
  > R.oo::attachLocally(dataS)
  > 
  > # Non-Paired PSCBS segmentation
  > fit <- segmentByNonPairedPSCBS(CT, betaT=betaT,
  +                             chromosome=chromosome, x=x,
  +                             seed=0xBEEF, verbose=verbose)
  > print(fit)
    chromosome tcnId dhId     start       end tcnNbrOfLoci tcnMean tcnNbrOfSNPs
  1          1     1    1    554484 143663981         1880  1.3916          778
  2          1     2    1 143663981 185240536          671  2.0925          275
  3          1     3    1 185240536 246679946         1111  2.6545          417
    tcnNbrOfHets dhNbrOfLoci    dhMean    c1Mean    c2Mean
  1          778         778 0.4009957 0.4167872 0.9748128
  2          275         275 0.2344486 0.8009582 1.2915418
  3          417         417 0.2819897 0.9529792 1.7015208
  > 
  > 
  > # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  > # Bootstrap segment level estimates
  > # (used by the AB caller, which, if skipped here,
  > #  will do it automatically)
  > # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  > fit <- bootstrapTCNandDHByRegion(fit, B=100, verbose=verbose)
  > print(fit)
    chromosome tcnId dhId     start       end tcnNbrOfLoci tcnMean tcnNbrOfSNPs
  1          1     1    1    554484 143663981         1880  1.3916          778
  2          1     2    1 143663981 185240536          671  2.0925          275
  3          1     3    1 185240536 246679946         1111  2.6545          417
    tcnNbrOfHets dhNbrOfLoci    dhMean    c1Mean    c2Mean
  1          778         778 0.4009957 0.4167872 0.9748128
  2          275         275 0.2344486 0.8009582 1.2915418
  3          417         417 0.2819897 0.9529792 1.7015208
  > 
  > 
  > # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  > # Calling segments in allelic balance (AB)
  > # NOTE: Ideally, this should be done on whole-genome data
  > # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  > # Explicitly estimate the threshold in DH for calling AB
  > # (which be done by default by the caller, if skipped here)
  > deltaAB <- estimateDeltaAB(fit, flavor="qq(DH)", verbose=verbose)
  Error in loadNamespace(name) : there is no package called ‘Hmisc’
  Calls: estimateDeltaAB ... loadNamespace -> withRestarts -> withOneRestart -> doWithOneRestart
  Execution halted

❯ checking tests ...
  See below...

❯ checking running R code from vignettes ...
    ‘CBS.tex.rsp’... OK
    ‘PairedPSCBS.tex.rsp’... failed
   WARNING
  Errors in running code in vignettes:
  when running code in ‘PairedPSCBS.tex.rsp’
    ...
   logical       9       1       1 
  Calling ROH...done
  > fit <- callROH(fit, verbose = -10)
  > withCapture({
  +     fit <- callAB(fit, verbose = -10)
  + })
  
    When sourcing ‘PairedPSCBS.R’:
  Error: there is no package called ‘Hmisc’
  Execution halted

❯ checking package dependencies ... NOTE
  Package suggested but not available for checking: ‘Hmisc’

❯ checking Rd cross-references ... NOTE
  Package unavailable to check Rd xrefs: ‘Hmisc’

❯ checking re-building of vignette outputs ... NOTE
  Error(s) in re-building vignettes:
  --- re-building ‘CBS.tex.rsp’ using rsp
  Segmenting by CBS...
   Chromosome: 1
   Segmenting multiple segments on current chromosome...
    Number of segments: 3
    Random seed temporarily set (seed=c(48879), kind="L'Ecuyer-CMRG")
    Produced 3 seeds from this stream for future usage
     Segmenting by CBS...
      Chromosome: 1
       Random seed temporarily set (seed=c(10407, 1066287653, -51199871, 161854402, -1995183193, 1503453565, -747102133), kind="L'Ecuyer-CMRG")
     Segmenting by CBS...done
     Segmenting by CBS...
      Chromosome: 1
       Random seed temporarily set (seed=c(10407, -821273412, -52578226, 1415511586, 721384351, -665928286, 1316562960), kind="L'Ecuyer-CMRG")
     Segmenting by CBS...done
   Segmenting multiple segments on current chromosome...done
  Segmenting by CBS...done
  Warning in plotTracks.CBS(fit) :
    Setting default 'Clim' assuming the signal type is ‘ratio’ because signalType(fit) is unknown (‘NA’). Use signalType(fit) <- ‘ratio’ to avoid this warning.
  Prune segments by hierarchical clustering...
   Clustering arguments:
   List of 4
    $ size        : NULL
    $ distMethod  : chr "euclidean"
    $ hclustMethod: chr "ward.D"
    $ h           : num 0.25
   Clustering...
    Hierarchical clustering of segmented copy numbers...
     Extracting/sampling CNs...
       num [1:9, 1] 1.38 3.19 1.39 NA 2.07 ...
       - attr(*, "dimnames")=List of 2
        ..$ : chr [1:9] "1" "2" "3" "4" ...
        ..$ : chr "mean"
       num [1:8, 1] 1.38 3.19 1.39 2.07 2.71 ...
       - attr(*, "dimnames")=List of 2
        ..$ : chr [1:8] "1" "2" "3" "5" ...
        ..$ : chr "mean"
     Extracting/sampling CNs...done
     Calculating distance matrix...
       'dist' num [1:28] 1.8055 0.0115 0.6858 1.3278 1.2065 ...
       - attr(*, "Size")= int 8
       - attr(*, "Labels")= chr [1:8] "1" "2" "3" "5" ...
       - attr(*, "Diag")= logi FALSE
       - attr(*, "Upper")= logi FALSE
       - attr(*, "method")= chr "euclidean"
       - attr(*, "call")= language stats::dist(x = C, method = distMethod)
     Calculating distance matrix...done
     Clustering...
      List of 7
       $ merge      : int [1:7, 1:2] -1 -6 -5 -2 -4 4 1 -3 -8 2 ...
       $ height     : num [1:7] 0.0115 0.0514 0.1103 0.6868 0.827 ...
       $ order      : int [1:8] 1 3 2 7 4 5 6 8
       $ labels     : chr [1:8] "1" "2" "3" "5" ...
       $ method     : chr "ward.D"
       $ call       : language stats::hclust(d = D, method = hclustMethod)
       $ dist.method: chr "euclidean"
       - attr(*, "class")= chr "hclust"
     Clustering...done
    Hierarchical clustering of segmented copy numbers...done
    
    Call:
    stats::hclust(d = D, method = hclustMethod)
    
    Cluster method   : ward.D 
    Distance         : euclidean 
    Number of objects: 8 
    
   Clustering...done
   Cutting tree...
    Cutting arguments:
    List of 2
     $ tree:List of 7
      ..$ merge      : int [1:7, 1:2] -1 -6 -5 -2 -4 4 1 -3 -8 2 ...
      ..$ height     : num [1:7] 0.0115 0.0514 0.1103 0.6868 0.827 ...
      ..$ order      : int [1:8] 1 3 2 7 4 5 6 8
      ..$ labels     : chr [1:8] "1" "2" "3" "5" ...
      ..$ method     : chr "ward.D"
      ..$ call       : language stats::hclust(d = D, method = hclustMethod)
      ..$ dist.method: chr "euclidean"
      ..- attr(*, "class")= chr "hclust"
     $ h   : num 0.25
     Named int [1:8] 1 2 1 3 4 4 5 4
     - attr(*, "names")= chr [1:8] "1" "2" "3" "5" ...
    List of 5
     $ 1: int [1:2] 1 3
     $ 2: int 2
     $ 3: int 5
     $ 4: int [1:3] 6 7 9
     $ 5: int 8
     - attr(*, "call")= language by.default(data = names(p), INDICES = p, FUN = functio..
     - attr(*, "class")= chr "by"
   Cutting tree...done
   Merging mean levels of clustered segments...
   Merging mean levels of clustered segments...done
   Merging neighboring segments within each cluster...
    Cluster #1 of 5...
     Segments in cluster:
      int [1:2] 1 3
     Left indices of neighboring segments:
      int(0) 
    Cluster #1 of 5...done
    Cluster #2 of 5...
     Segments in cluster:
      int 2
     Left indices of neighboring segments:
      int(0) 
    Cluster #2 of 5...done
    Cluster #3 of 5...
     Segments in cluster:
      int 5
     Left indices of neighboring segments:
      int(0) 
    Cluster #3 of 5...done
    Cluster #4 of 5...
     Segments in cluster:
      int [1:3] 6 7 9
     Left indices of neighboring segments:
      int 6
    Cluster #4 of 5...done
    Cluster #5 of 5...
     Segments in cluster:
      int 8
     Left indices of neighboring segments:
      int(0) 
    Cluster #5 of 5...done
    Left indices of segments to be merged:
     int 6
   Merging neighboring segments within each cluster...done
   Merging segments...
   Merging segments...done
   Updating segment means...
   Updating segment means...done
  Prune segments by hierarchical clustering...done
  Warning in plotTracks.CBS(fitP) :
    Setting default 'Clim' assuming the signal type is ‘ratio’ because signalType(fitP) is unknown (‘NA’). Use signalType(fitP) <- ‘ratio’ to avoid this warning.
  --- finished re-building ‘CBS.tex.rsp’
  
  --- re-building ‘PairedPSCBS.tex.rsp’ using rsp
  Segmenting paired tumor-normal signals using Paired PSCBS...
   Calling genotypes from normal allele B fractions...
     num [1:73346] 0.827 0.875 0.887 0.884 0.103 ...
    Called genotypes:
     num [1:73346] 1 1 1 1 0 0.5 0.5 0 0.5 0 ...
     - attr(*, "modelFit")=List of 1
      ..$ :List of 7
      .. ..$ flavor             : chr "density"
      .. ..$ cn                 : int 2
      .. ..$ nbrOfGenotypeGroups: int 3
      .. ..$ tau                : num [1:2] 0.316 0.679
      .. ..$ n                  : int 73181
      .. ..$ fit                :Classes ‘PeaksAndValleys’ and 'data.frame':	5 obs. of ..
      .. .. ..$ type   : chr [1:5] "peak" "valley" "peak" "valley" ...
      .. .. ..$ x      : num [1:5] 0.096 0.316 0.496 0.679 0.894
      .. .. ..$ density: num [1:5] 1.709 0.437 1.12 0.466 1.673
      .. ..$ fitValleys         :Classes ‘PeaksAndValleys’ and 'data.frame':	2 obs. of ..
      .. .. ..$ type   : chr [1:2] "valley" "valley"
      .. .. ..$ x      : num [1:2] 0.316 0.679
      .. .. ..$ density: num [1:2] 0.437 0.466
      ..- attr(*, "class")= chr [1:2] "NaiveGenotypeModelFit" "list"
    muN
        0   0.5     1 
    26499 20487 26360 
   Calling genotypes from normal allele B fractions...done
   Normalizing betaT using betaN (TumorBoost)...
    Normalized BAFs:
     num [1:73346] 0.93 0.896 0.947 0.894 0.126 ...
     - attr(*, "modelFit")=List of 5
      ..$ method       : chr "normalizeTumorBoost"
      ..$ flavor       : chr "v4"
      ..$ delta        : num [1:73346] -0.173 -0.125 -0.113 -0.116 0.103 ...
      .. ..- attr(*, "modelFit")=List of 1
      .. .. ..$ :List of 7
      .. .. .. ..$ flavor             : chr "density"
      .. .. .. ..$ cn                 : int 2
      .. .. .. ..$ nbrOfGenotypeGroups: int 3
      .. .. .. ..$ tau                : num [1:2] 0.316 0.679
      .. .. .. ..$ n                  : int 73181
      .. .. .. ..$ fit                :Classes ‘PeaksAndValleys’ and 'data.frame':	5 ob..
      .. .. .. .. ..$ type   : chr [1:5] "peak" "valley" "peak" "valley" ...
      .. .. .. .. ..$ x      : num [1:5] 0.096 0.316 0.496 0.679 0.894
      .. .. .. .. ..$ density: num [1:5] 1.709 0.437 1.12 0.466 1.673
      .. .. .. ..$ fitValleys         :Classes ‘PeaksAndValleys’ and 'data.frame':	2 ob..
      .. .. .. .. ..$ type   : chr [1:2] "valley" "valley"
      .. .. .. .. ..$ x      : num [1:2] 0.316 0.679
      .. .. .. .. ..$ density: num [1:2] 0.437 0.466
      .. .. ..- attr(*, "class")= chr [1:2] "NaiveGenotypeModelFit" "list"
      ..$ preserveScale: logi FALSE
      ..$ scaleFactor  : num NA
   Normalizing betaT using betaN (TumorBoost)...done
   Setup up data...
    'data.frame':	73346 obs. of  7 variables:
     $ chromosome: int  1 1 1 1 1 1 1 1 1 1 ...
     $ x         : num  1145994 2224111 2319424 2543484 2926730 ...
     $ CT        : num  1.625 1.071 1.406 1.18 0.856 ...
     $ betaT     : num  0.757 0.771 0.834 0.778 0.229 ...
     $ betaTN    : num  0.93 0.896 0.947 0.894 0.126 ...
      ..- attr(*, "modelFit")=List of 5
      .. ..$ method       : chr "normalizeTumorBoost"
      .. ..$ flavor       : chr "v4"
      .. ..$ delta        : num [1:73346] -0.173 -0.125 -0.113 -0.116 0.103 ...
      .. .. ..- attr(*, "modelFit")=List of 1
      .. .. .. ..$ :List of 7
      .. .. .. .. ..$ flavor             : chr "density"
      .. .. .. .. ..$ cn                 : int 2
      .. .. .. .. ..$ nbrOfGenotypeGroups: int 3
      .. .. .. .. ..$ tau                : num [1:2] 0.316 0.679
      .. .. .. .. ..$ n                  : int 73181
      .. .. .. .. ..$ fit                :Classes ‘PeaksAndValleys’ and 'data.frame':	5..
      .. .. .. .. .. ..$ type   : chr [1:5] "peak" "valley" "peak" "valley" ...
      .. .. .. .. .. ..$ x      : num [1:5] 0.096 0.316 0.496 0.679 0.894
      .. .. .. .. .. ..$ density: num [1:5] 1.709 0.437 1.12 0.466 1.673
      .. .. .. .. ..$ fitValleys         :Classes ‘PeaksAndValleys’ and 'data.frame':	2..
      .. .. .. .. .. ..$ type   : chr [1:2] "valley" "valley"
      .. .. .. .. .. ..$ x      : num [1:2] 0.316 0.679
      .. .. .. .. .. ..$ density: num [1:2] 0.437 0.466
      .. .. .. ..- attr(*, "class")= chr [1:2] "NaiveGenotypeModelFit" "list"
      .. ..$ preserveScale: logi FALSE
      .. ..$ scaleFactor  : num NA
     $ betaN     : num  0.827 0.875 0.887 0.884 0.103 ...
     $ muN       : num  1 1 1 1 0 0.5 0.5 0 0.5 0 ...
      ..- attr(*, "modelFit")=List of 1
      .. ..$ :List of 7
      .. .. ..$ flavor             : chr "density"
      .. .. ..$ cn                 : int 2
      .. .. ..$ nbrOfGenotypeGroups: int 3
      .. .. ..$ tau                : num [1:2] 0.316 0.679
      .. .. ..$ n                  : int 73181
      .. .. ..$ fit                :Classes ‘PeaksAndValleys’ and 'data.frame':	5 obs. ..
      .. .. .. ..$ type   : chr [1:5] "peak" "valley" "peak" "valley" ...
      .. .. .. ..$ x      : num [1:5] 0.096 0.316 0.496 0.679 0.894
      .. .. .. ..$ density: num [1:5] 1.709 0.437 1.12 0.466 1.673
      .. .. ..$ fitValleys         :Classes ‘PeaksAndValleys’ and 'data.frame':	2 obs. ..
      .. .. .. ..$ type   : chr [1:2] "valley" "valley"
      .. .. .. ..$ x      : num [1:2] 0.316 0.679
      .. .. .. ..$ density: num [1:2] 0.437 0.466
      .. ..- attr(*, "class")= chr [1:2] "NaiveGenotypeModelFit" "list"
   Setup up data...done
   Dropping loci for which TCNs are missing...
    Number of loci dropped: 49
   Dropping loci for which TCNs are missing...done
   Ordering data along genome...
    'data.frame':	73297 obs. of  7 variables:
     $ chromosome: int  1 1 1 1 1 1 1 1 1 1 ...
     $ x         : num  554484 554636 557616 711153 730720 ...
     $ CT        : num  1.878 0.166 2.119 1.397 1.799 ...
     $ betaT     : num  0.0646 0.6462 0.8257 0.2598 0.1672 ...
     $ betaTN    : num  -0.0515 0.9543 0.9831 0.0832 -0.1172 ...
     $ betaN     : num  0.116 0.692 0.843 0.177 0.284 ...
     $ muN       : num  0 1 1 0 0 1 0.5 0 1 1 ...
   Ordering data along genome...done
   Keeping only current chromosome for 'knownSegments'...
    Chromosome: 1
    Known segments for this chromosome:
      chromosome    start      end   length
    1          1     -Inf 1.21e+08      Inf
    2          1 1.21e+08 1.42e+08 20517398
    3          1 1.42e+08      Inf      Inf
   Keeping only current chromosome for 'knownSegments'...done
   alphaTCN: 0.009
   alphaDH: 0.001
   Number of loci: 73297
   Calculating DHs...
    Number of SNPs: 73297
    Number of heterozygous SNPs: 20472 (27.93%)
    Normalized DHs:
     num [1:73297] NA NA NA NA NA ...
   Calculating DHs...done
   Random seed temporarily set (seed=c(48879), kind="L'Ecuyer-CMRG")
   Produced 2 seeds from this stream for future usage
   Identification of change points by total copy numbers...
    Segmenting by CBS...
     Chromosome: 1
     Segmenting multiple segments on current chromosome...
      Number of segments: 3
      Random seed temporarily set (seed=c(10407, 1066287653, -51199871, 161854402, -1995183193, 1503453565, -747102133), kind="L'Ecuyer-CMRG")
      Produced 3 seeds from this stream for future usage
       Segmenting by CBS...
        Chromosome: 1
         Random seed temporarily set (seed=c(10407, 1797822437, 388243314, 91406689, -519578635, -1381924756, 1089253019), kind="L'Ecuyer-CMRG")
       Segmenting by CBS...done
       Segmenting by CBS...
        Chromosome: 1
         Random seed temporarily set (seed=c(10407, -1924040949, -1632234809, -437763632, -1464377300, 676654412, 2080370711), kind="L'Ecuyer-CMRG")
       Segmenting by CBS...done
     Segmenting multiple segments on current chromosome...done
    Segmenting by CBS...done
    List of 4
     $ data   :'data.frame':	73297 obs. of  4 variables:
      ..$ chromosome: int [1:73297] 1 1 1 1 1 1 1 1 1 1 ...
      ..$ x         : num [1:73297] 554484 554636 557616 711153 730720 ...
      ..$ y         : num [1:73297] 1.878 0.166 2.119 1.397 1.799 ...
      ..$ index     : int [1:73297] 1 2 3 4 5 6 7 8 9 10 ...
     $ output :'data.frame':	7 obs. of  6 variables:
      ..$ sampleName: chr [1:7] NA NA NA NA ...
      ..$ chromosome: int [1:7] 1 1 1 1 1 1 1
      ..$ start     : num [1:7] 5.54e+05 1.21e+08 1.42e+08 1.86e+08 1.99e+08 ...
      ..$ end       : num [1:7] 1.21e+08 1.42e+08 1.86e+08 1.99e+08 2.07e+08 ...
      ..$ nbrOfLoci : int [1:7] 37495 0 13434 4018 2755 14 15581
      ..$ mean      : num [1:7] 1.38 NA 2.07 2.71 2.59 ...
     $ segRows:'data.frame':	7 obs. of  2 variables:
      ..$ startRow: int [1:7] 1 NA 37496 50930 54948 57703 57717
      ..$ endRow  : int [1:7] 37495 NA 50929 54947 57702 57716 73297
     $ params :List of 5
      ..$ alpha        : num 0.009
      ..$ undo         : num 0
      ..$ joinSegments : logi TRUE
      ..$ knownSegments:'data.frame':	4 obs. of  3 variables:
      .. ..$ chromosome: int [1:4] 1 1 2 1
      .. ..$ start     : num [1:4] -Inf -Inf -Inf 1.42e+08
      .. ..$ end       : num [1:4] 1.21e+08 Inf Inf Inf
      ..$ seed         : int [1:7] 10407 1066287653 -51199871 161854402 -1995183193 150..
     - attr(*, "class")= chr [1:2] "CBS" "AbstractCBS"
     - attr(*, "processingTime")= 'proc_time' Named num [1:5] 8.78 0 8.78 0 0
      ..- attr(*, "names")= chr [1:5] "user.self" "sys.self" "elapsed" "user.child" ...
     - attr(*, "pkgDetails")= chr "DNAcopy v1.64.0"
     - attr(*, "randomSeed")= int [1:7] 10407 1797822437 388243314 91406689 -519578635 ..
   Identification of change points by total copy numbers...done
   Restructure TCN segmentation results...
      chromosome tcnStart   tcnEnd tcnNbrOfLoci tcnMean
    1          1 5.54e+05 1.21e+08        37495    1.38
    2          1 1.21e+08 1.42e+08            0      NA
    3          1 1.42e+08 1.86e+08        13434    2.07
    4          1 1.86e+08 1.99e+08         4018    2.71
    5          1 1.99e+08 2.07e+08         2755    2.59
    6          1 2.07e+08 2.07e+08           14    3.87
    7          1 2.07e+08 2.47e+08        15581    2.64
    Number of TCN segments: 7
   Restructure TCN segmentation results...done
   Total CN segment #1 ([    554484,1.20993e+08]) of 7...
    Number of TCN loci in segment: 37495
    Locus data for TCN segment:
    'data.frame':	37495 obs. of  9 variables:
     $ chromosome: int  1 1 1 1 1 1 1 1 1 1 ...
     $ x         : num  554484 554636 557616 711153 730720 ...
     $ CT        : num  1.878 0.166 2.119 1.397 1.799 ...
     $ betaT     : num  0.0646 0.6462 0.8257 0.2598 0.1672 ...
     $ betaTN    : num  -0.0515 0.9543 0.9831 0.0832 -0.1172 ...
     $ betaN     : num  0.116 0.692 0.843 0.177 0.284 ...
     $ muN       : num  0 1 1 0 0 1 0.5 0 1 1 ...
     $ index     : int  1 2 3 4 5 6 7 8 9 10 ...
     $ rho       : num  NA NA NA NA NA ...
    Number of loci: 37495
    Number of SNPs: 10146 (27.06%)
    Number of heterozygous SNPs: 10146 (100.00%)
    Chromosome: 1
    Segmenting DH signals...
     Segmenting by CBS...
      Chromosome: 1
       Random seed temporarily set (seed=c(10407, 1797822437, 388243314, 91406689, -519578635, -1381924756, 1089253019), kind="L'Ecuyer-CMRG")
     Segmenting by CBS...done
     List of 4
      $ data   :'data.frame':	37495 obs. of  4 variables:
       ..$ chromosome: int [1:37495] 1 1 1 1 1 1 1 1 1 1 ...
       ..$ x         : num [1:37495] 554484 554636 557616 711153 730720 ...
       ..$ y         : num [1:37495] NA NA NA NA NA ...
       ..$ index     : int [1:37495] 1 2 3 4 5 6 7 8 9 10 ...
      $ output :'data.frame':	5 obs. of  6 variables:
       ..$ sampleName: chr [1:5] NA NA NA NA ...
       ..$ chromosome: int [1:5] 1 1 1 1 1
       ..$ start     : num [1:5] 5.54e+05 4.28e+06 4.27e+07 1.20e+08 1.20e+08
       ..$ end       : num [1:5] 4.28e+06 4.27e+07 1.20e+08 1.20e+08 1.21e+08
       ..$ nbrOfLoci : int [1:5] 169 3260 6657 8 52
       ..$ mean      : num [1:5] 0.4145 0.4944 0.5205 0.0767 0.5123
      $ segRows:'data.frame':	5 obs. of  2 variables:
       ..$ startRow: int [1:5] 7 693 12094 37253 37332
       ..$ endRow  : int [1:5] 684 12093 37250 37320 37449
      $ params :List of 5
       ..$ alpha        : num 0.001
       ..$ undo         : num 0
       ..$ joinSegments : logi TRUE
       ..$ knownSegments:'data.frame':	1 obs. of  3 variables:
       .. ..$ chromosome: int 1
       .. ..$ start     : num 554484
       .. ..$ end       : num 1.21e+08
       ..$ seed         : int [1:7] 10407 1797822437 388243314 91406689 -519578635 -1381..
      - attr(*, "class")= chr [1:2] "CBS" "AbstractCBS"
      - attr(*, "processingTime")= 'proc_time' Named num [1:5] 1.91 0.004 1.914 0 0
       ..- attr(*, "names")= chr [1:5] "user.self" "sys.self" "elapsed" "user.child" ...
      - attr(*, "pkgDetails")= chr "DNAcopy v1.64.0"
      - attr(*, "randomSeed")= int [1:7] 10407 1797822437 388243314 91406689 -519578635 ..
     DH segmentation (locally-indexed) rows:
       startRow endRow
     1        7    684
     2      693  12093
     3    12094  37250
     4    37253  37320
     5    37332  37449
      int [1:37495] 1 2 3 4 5 6 7 8 9 10 ...
     DH segmentation rows:
       startRow endRow
     1        7    684
     2      693  12093
     3    12094  37250
     4    37253  37320
     5    37332  37449
    Segmenting DH signals...done
    DH segmentation table:
       dhStart    dhEnd dhNbrOfLoci dhMean
    1 5.54e+05 4.28e+06         169 0.4145
    2 4.28e+06 4.27e+07        3260 0.4944
    3 4.27e+07 1.20e+08        6657 0.5205
    4 1.20e+08 1.20e+08           8 0.0767
    5 1.20e+08 1.21e+08          52 0.5123
      startRow endRow
    1        7    684
    2      693  12093
    3    12094  37250
    4    37253  37320
    5    37332  37449
    Rows:
    [1] 1 1 1 1 1
    TCN segmentation rows:
        startRow endRow
    1          1  37495
    1.1        1  37495
    1.2        1  37495
    1.3        1  37495
    1.4        1  37495
    TCN and DH segmentation rows:
        startRow endRow
    1          1  37495
    1.1        1  37495
    1.2        1  37495
    1.3        1  37495
    1.4        1  37495
      startRow endRow
    1        7    684
    2      693  12093
    3    12094  37250
    4    37253  37320
    5    37332  37449
    NULL
    TCN segmentation (expanded) rows:
        startRow endRow
    1          1  37495
    1.1        1  37495
    1.2        1  37495
    1.3        1  37495
    1.4        1  37495
    TCN and DH segmentation rows:
      startRow endRow
    1        1  37495
    2       NA     NA
    3    37496  50929
    4    50930  54947
    5    54948  57702
    6    57703  57716
    7    57717  73297
      startRow endRow
    1        7    684
    2      693  12093
    3    12094  37250
    4    37253  37320
    5    37332  37449
      startRow endRow
    1        1  37495
    2        1  37495
    3        1  37495
    4        1  37495
    5        1  37495
    Total CN segmentation table (expanded):
        chromosome tcnStart   tcnEnd tcnNbrOfLoci tcnMean tcnNbrOfSNPs tcnNbrOfHets
    1            1   554484 1.21e+08        37495    1.38        10146        10146
    1.1          1   554484 1.21e+08        37495    1.38        10146        10146
    1.2          1   554484 1.21e+08        37495    1.38        10146        10146
    1.3          1   554484 1.21e+08        37495    1.38        10146        10146
    1.4          1   554484 1.21e+08        37495    1.38        10146        10146
    (TCN,DH) segmentation for one total CN segment:
        tcnId dhId chromosome tcnStart   tcnEnd tcnNbrOfLoci tcnMean tcnNbrOfSNPs
    1       1    1          1   554484 1.21e+08        37495    1.38        10146
    1.1     1    2          1   554484 1.21e+08        37495    1.38        10146
    1.2     1    3          1   554484 1.21e+08        37495    1.38        10146
    1.3     1    4          1   554484 1.21e+08        37495    1.38        10146
    1.4     1    5          1   554484 1.21e+08        37495    1.38        10146
        tcnNbrOfHets  dhStart    dhEnd dhNbrOfLoci dhMean
    1          10146 5.54e+05 4.28e+06         169 0.4145
    1.1        10146 4.28e+06 4.27e+07        3260 0.4944
    1.2        10146 4.27e+07 1.20e+08        6657 0.5205
    1.3        10146 1.20e+08 1.20e+08           8 0.0767
    1.4        10146 1.20e+08 1.21e+08          52 0.5123
   Total CN segment #1 ([    554484,1.20993e+08]) of 7...done
   Total CN segment #2 ([1.20993e+08,1.4151e+08]) of 7...
    Number of TCN loci in segment: 0
    Locus data for TCN segment:
    'data.frame':	0 obs. of  9 variables:
     $ chromosome: int 
     $ x         : num 
     $ CT        : num 
     $ betaT     : num 
     $ betaTN    : num 
     $ betaN     : num 
     $ muN       : num 
     $ index     : int 
     $ rho       : num 
    Number of loci: 0
    Number of SNPs: 0 (NaN%)
    Number of heterozygous SNPs: 0 (NaN%)
    Chromosome: 1
    Segmenting DH signals...
     Segmenting by CBS...
      Chromosome: NA
       Random seed temporarily set (seed=c(10407, 1797822437, 388243314, 91406689, -519578635, -1381924756, 1089253019), kind="L'Ecuyer-CMRG")
     Segmenting by CBS...done
     List of 4
      $ data   :'data.frame':	0 obs. of  4 variables:
       ..$ chromosome: int(0) 
       ..$ x         : num(0) 
       ..$ y         : num(0) 
       ..$ index     : int(0) 
      $ output :'data.frame':	0 obs. of  6 variables:
       ..$ sampleName: chr(0) 
       ..$ chromosome: num(0) 
       ..$ start     : num(0) 
       ..$ end       : num(0) 
       ..$ nbrOfLoci : int(0) 
       ..$ mean      : num(0) 
      $ segRows:'data.frame':	0 obs. of  2 variables:
       ..$ startRow: int(0) 
       ..$ endRow  : int(0) 
      $ params :List of 5
       ..$ alpha        : num 0.001
       ..$ undo         : num 0
       ..$ joinSegments : logi TRUE
       ..$ knownSegments:'data.frame':	0 obs. of  3 variables:
       .. ..$ chromosome: int(0) 
       .. ..$ start     : num(0) 
       .. ..$ end       : num(0) 
       ..$ seed         : int [1:7] 10407 1797822437 388243314 91406689 -519578635 -1381..
      - attr(*, "class")= chr [1:2] "CBS" "AbstractCBS"
      - attr(*, "processingTime")= 'proc_time' Named num [1:5] 0.002 0 0.002 0 0
       ..- attr(*, "names")= chr [1:5] "user.self" "sys.self" "elapsed" "user.child" ...
      - attr(*, "pkgDetails")= chr "DNAcopy v1.64.0"
      - attr(*, "randomSeed")= int [1:7] 10407 1797822437 388243314 91406689 -519578635 ..
     DH segmentation (locally-indexed) rows:
     [1] startRow endRow  
     <0 rows> (or 0-length row.names)
      int(0) 
     DH segmentation rows:
     [1] startRow endRow  
     <0 rows> (or 0-length row.names)
    Segmenting DH signals...done
    DH segmentation table:
       dhStart dhEnd dhNbrOfLoci dhMean
    NA      NA    NA          NA     NA
       startRow endRow
    NA       NA     NA
    Rows:
    [1] 2
    TCN segmentation rows:
      startRow endRow
    2       NA     NA
    TCN and DH segmentation rows:
      startRow endRow
    2       NA     NA
       startRow endRow
    NA       NA     NA
      startRow endRow
    1        1  37495
    2        1  37495
    3        1  37495
    4        1  37495
    5        1  37495
    TCN segmentation (expanded) rows:
       startRow endRow
    1         1  37495
    2         1  37495
    3         1  37495
    4         1  37495
    5         1  37495
    21       NA     NA
    TCN and DH segmentation rows:
      startRow endRow
    1        1  37495
    2       NA     NA
    3    37496  50929
    4    50930  54947
    5    54948  57702
    6    57703  57716
    7    57717  73297
      startRow endRow
    1        7    684
    2      693  12093
    3    12094  37250
    4    37253  37320
    5    37332  37449
    6       NA     NA
      startRow endRow
    1        1  37495
    2        1  37495
    3        1  37495
    4        1  37495
    5        1  37495
    6       NA     NA
    Total CN segmentation table (expanded):
      chromosome tcnStart   tcnEnd tcnNbrOfLoci tcnMean tcnNbrOfSNPs tcnNbrOfHets
    2          1 1.21e+08 1.42e+08            0      NA            0            0
    (TCN,DH) segmentation for one total CN segment:
      tcnId dhId chromosome tcnStart   tcnEnd tcnNbrOfLoci tcnMean tcnNbrOfSNPs
    2     2    1          1 1.21e+08 1.42e+08            0      NA            0
      tcnNbrOfHets dhStart dhEnd dhNbrOfLoci dhMean
    2            0      NA    NA          NA     NA
   Total CN segment #2 ([1.20993e+08,1.4151e+08]) of 7...done
   Total CN segment #3 ([1.4151e+08,1.85528e+08]) of 7...
    Number of TCN loci in segment: 13434
    Locus data for TCN segment:
    'data.frame':	13434 obs. of  9 variables:
     $ chromosome: int  1 1 1 1 1 1 1 1 1 1 ...
     $ x         : num  1.42e+08 1.42e+08 1.42e+08 1.42e+08 1.42e+08 ...
     $ CT        : num  1.74 1.98 3.48 2.11 2.48 ...
     $ betaT     : num  0.0282 0.9272 0.9103 0.0544 0.5531 ...
     $ betaTN    : num  -0.02304 0.96666 0.97572 0.00855 0.48389 ...
     $ betaN     : num  0.0512 0.9605 0.9346 0.0459 0.5715 ...
     $ muN       : num  0 1 1 0 0.5 0.5 0 0.5 0.5 1 ...
     $ index     : int  37496 37497 37498 37499 37500 37501 37502 37503 37504 37505 ...
     $ rho       : num  NA NA NA NA 0.0322 ...
    Number of loci: 13434
    Number of SNPs: 3770 (28.06%)
    Number of heterozygous SNPs: 3770 (100.00%)
    Chromosome: 1
    Segmenting DH signals...
     Segmenting by CBS...
      Chromosome: 1
       Random seed temporarily set (seed=c(10407, 1797822437, 388243314, 91406689, -519578635, -1381924756, 1089253019), kind="L'Ecuyer-CMRG")
     Segmenting by CBS...done
     List of 4
      $ data   :'data.frame':	13434 obs. of  4 variables:
       ..$ chromosome: int [1:13434] 1 1 1 1 1 1 1 1 1 1 ...
       ..$ x         : num [1:13434] 1.42e+08 1.42e+08 1.42e+08 1.42e+08 1.42e+08 ...
       ..$ y         : num [1:13434] NA NA NA NA 0.0322 ...
       ..$ index     : int [1:13434] 1 2 3 4 5 6 7 8 9 10 ...
      $ output :'data.frame':	1 obs. of  6 variables:
       ..$ sampleName: chr NA
       ..$ chromosome: int 1
       ..$ start     : num 1.42e+08
       ..$ end       : num 1.86e+08
       ..$ nbrOfLoci : int 3770
       ..$ mean      : num 0.0943
      $ segRows:'data.frame':	1 obs. of  2 variables:
       ..$ startRow: int 5
       ..$ endRow  : int 13434
      $ params :List of 5
       ..$ alpha        : num 0.001
       ..$ undo         : num 0
       ..$ joinSegments : logi TRUE
       ..$ knownSegments:'data.frame':	1 obs. of  3 variables:
       .. ..$ chromosome: int 1
       .. ..$ start     : num 1.42e+08
       .. ..$ end       : num 1.86e+08
       ..$ seed         : int [1:7] 10407 1797822437 388243314 91406689 -519578635 -1381..
      - attr(*, "class")= chr [1:2] "CBS" "AbstractCBS"
      - attr(*, "processingTime")= 'proc_time' Named num [1:5] 0.102 0 0.102 0 0
       ..- attr(*, "names")= chr [1:5] "user.self" "sys.self" "elapsed" "user.child" ...
      - attr(*, "pkgDetails")= chr "DNAcopy v1.64.0"
      - attr(*, "randomSeed")= int [1:7] 10407 1797822437 388243314 91406689 -519578635 ..
     DH segmentation (locally-indexed) rows:
       startRow endRow
     1        5  13434
      int [1:13434] 37496 37497 37498 37499 37500 37501 37502 37503 37504 37505 ...
     DH segmentation rows:
       startRow endRow
     1    37500  50929
    Segmenting DH signals...done
    DH segmentation table:
       dhStart    dhEnd dhNbrOfLoci dhMean
    1 1.42e+08 1.86e+08        3770 0.0943
      startRow endRow
    1    37500  50929
    Rows:
    [1] 3
    TCN segmentation rows:
      startRow endRow
    3    37496  50929
    TCN and DH segmentation rows:
      startRow endRow
    3    37496  50929
      startRow endRow
    1    37500  50929
      startRow endRow
    1        1  37495
    2        1  37495
    3        1  37495
    4        1  37495
    5        1  37495
    6       NA     NA
    TCN segmentation (expanded) rows:
       startRow endRow
    1         1  37495
    2         1  37495
    3         1  37495
    4         1  37495
    5         1  37495
    6        NA     NA
    31    37496  50929
    TCN and DH segmentation rows:
      startRow endRow
    1        1  37495
    2       NA     NA
    3    37496  50929
    4    50930  54947
    5    54948  57702
    6    57703  57716
    7    57717  73297
      startRow endRow
    1        7    684
    2      693  12093
    3    12094  37250
    4    37253  37320
    5    37332  37449
    6       NA     NA
    7    37500  50929
      startRow endRow
    1        1  37495
    2        1  37495
    3        1  37495
    4        1  37495
    5        1  37495
    6       NA     NA
    7    37496  50929
    Total CN segmentation table (expanded):
      chromosome tcnStart   tcnEnd tcnNbrOfLoci tcnMean tcnNbrOfSNPs tcnNbrOfHets
    3          1 1.42e+08 1.86e+08        13434    2.07         3770         3770
    (TCN,DH) segmentation for one total CN segment:
      tcnId dhId chromosome tcnStart   tcnEnd tcnNbrOfLoci tcnMean tcnNbrOfSNPs
    3     3    1          1 1.42e+08 1.86e+08        13434    2.07         3770
      tcnNbrOfHets  dhStart    dhEnd dhNbrOfLoci dhMean
    3         3770 1.42e+08 1.86e+08        3770 0.0943
   Total CN segment #3 ([1.4151e+08,1.85528e+08]) of 7...done
   Total CN segment #4 ([1.85528e+08,1.99122e+08]) of 7...
    Number of TCN loci in segment: 4018
    Locus data for TCN segment:
    'data.frame':	4018 obs. of  9 variables:
     $ chromosome: int  1 1 1 1 1 1 1 1 1 1 ...
     $ x         : num  1.86e+08 1.86e+08 1.86e+08 1.86e+08 1.86e+08 ...
     $ CT        : num  2.73 2.46 2.17 3.41 2.43 ...
     $ betaT     : num  0.499 0.751 0.147 0.954 0.197 ...
     $ betaTN    : num  0.5496 0.9681 -0.0136 1.042 0.2571 ...
     $ betaN     : num  0.444 0.783 0.16 0.912 0.383 ...
     $ muN       : num  0.5 1 0 1 0.5 1 1 1 0 1 ...
     $ index     : int  50930 50931 50932 50933 50934 50935 50936 50937 50938 50939 ...
     $ rho       : num  0.0992 NA NA NA 0.4857 ...
    Number of loci: 4018
    Number of SNPs: 1271 (31.63%)
    Number of heterozygous SNPs: 1271 (100.00%)
    Chromosome: 1
    Segmenting DH signals...
     Segmenting by CBS...
      Chromosome: 1
       Random seed temporarily set (seed=c(10407, 1797822437, 388243314, 91406689, -519578635, -1381924756, 1089253019), kind="L'Ecuyer-CMRG")
     Segmenting by CBS...done
     List of 4
      $ data   :'data.frame':	4018 obs. of  4 variables:
       ..$ chromosome: int [1:4018] 1 1 1 1 1 1 1 1 1 1 ...
       ..$ x         : num [1:4018] 1.86e+08 1.86e+08 1.86e+08 1.86e+08 1.86e+08 ...
       ..$ y         : num [1:4018] 0.0992 NA NA NA 0.4857 ...
       ..$ index     : int [1:4018] 1 2 3 4 5 6 7 8 9 10 ...
      $ output :'data.frame':	1 obs. of  6 variables:
       ..$ sampleName: chr NA
       ..$ chromosome: int 1
       ..$ start     : num 1.86e+08
       ..$ end       : num 1.99e+08
       ..$ nbrOfLoci : int 1271
       ..$ mean      : num 0.256
      $ segRows:'data.frame':	1 obs. of  2 variables:
       ..$ startRow: int 1
       ..$ endRow  : int 4011
      $ params :List of 5
       ..$ alpha        : num 0.001
       ..$ undo         : num 0
       ..$ joinSegments : logi TRUE
       ..$ knownSegments:'data.frame':	1 obs. of  3 variables:
       .. ..$ chromosome: int 1
       .. ..$ start     : num 1.86e+08
       .. ..$ end       : num 1.99e+08
       ..$ seed         : int [1:7] 10407 1797822437 388243314 91406689 -519578635 -1381..
      - attr(*, "class")= chr [1:2] "CBS" "AbstractCBS"
      - attr(*, "processingTime")= 'proc_time' Named num [1:5] 0.016 0 0.016 0 0
       ..- attr(*, "names")= chr [1:5] "user.self" "sys.self" "elapsed" "user.child" ...
      - attr(*, "pkgDetails")= chr "DNAcopy v1.64.0"
      - attr(*, "randomSeed")= int [1:7] 10407 1797822437 388243314 91406689 -519578635 ..
     DH segmentation (locally-indexed) rows:
       startRow endRow
     1        1   4011
      int [1:4018] 50930 50931 50932 50933 50934 50935 50936 50937 50938 50939 ...
     DH segmentation rows:
       startRow endRow
     1    50930  54940
    Segmenting DH signals...done
    DH segmentation table:
       dhStart    dhEnd dhNbrOfLoci dhMean
    1 1.86e+08 1.99e+08        1271  0.256
      startRow endRow
    1    50930  54940
    Rows:
    [1] 4
    TCN segmentation rows:
      startRow endRow
    4    50930  54947
    TCN and DH segmentation rows:
      startRow endRow
    4    50930  54947
      startRow endRow
    1    50930  54940
      startRow endRow
    1        1  37495
    2        1  37495
    3        1  37495
    4        1  37495
    5        1  37495
    6       NA     NA
    7    37496  50929
    TCN segmentation (expanded) rows:
       startRow endRow
    1         1  37495
    2         1  37495
    3         1  37495
    4         1  37495
    5         1  37495
    6        NA     NA
    7     37496  50929
    41    50930  54947
    TCN and DH segmentation rows:
      startRow endRow
    1        1  37495
    2       NA     NA
    3    37496  50929
    4    50930  54947
    5    54948  57702
    6    57703  57716
    7    57717  73297
      startRow endRow
    1        7    684
    2      693  12093
    3    12094  37250
    4    37253  37320
    5    37332  37449
    6       NA     NA
    7    37500  50929
    8    50930  54940
      startRow endRow
    1        1  37495
    2        1  37495
    3        1  37495
    4        1  37495
    5        1  37495
    6       NA     NA
    7    37496  50929
    8    50930  54947
    Total CN segmentation table (expanded):
      chromosome tcnStart   tcnEnd tcnNbrOfLoci tcnMean tcnNbrOfSNPs tcnNbrOfHets
    4          1 1.86e+08 1.99e+08         4018    2.71         1271         1271
    (TCN,DH) segmentation for one total CN segment:
      tcnId dhId chromosome tcnStart   tcnEnd tcnNbrOfLoci tcnMean tcnNbrOfSNPs
    4     4    1          1 1.86e+08 1.99e+08         4018    2.71         1271
      tcnNbrOfHets  dhStart    dhEnd dhNbrOfLoci dhMean
    4         1271 1.86e+08 1.99e+08        1271  0.256
   Total CN segment #4 ([1.85528e+08,1.99122e+08]) of 7...done
   Total CN segment #5 ([1.99122e+08,2.06513e+08]) of 7...
    Number of TCN loci in segment: 2755
    Locus data for TCN segment:
    'data.frame':	2755 obs. of  9 variables:
     $ chromosome: int  1 1 1 1 1 1 1 1 1 1 ...
     $ x         : num  1.99e+08 1.99e+08 1.99e+08 1.99e+08 1.99e+08 ...
     $ CT        : num  2.15 2.75 2.18 2.27 2.52 ...
     $ betaT     : num  0.908 0.185 0.806 0.377 0.357 ...
     $ betaTN    : num  0.934 -0.033 0.956 0.401 0.354 ...
     $ betaN     : num  0.974 0.218 0.85 0.47 0.504 ...
     $ muN       : num  1 0 1 0.5 0.5 0.5 0.5 0.5 0 1 ...
     $ index     : int  54948 54949 54950 54951 54952 54953 54954 54955 54956 54957 ...
     $ rho       : num  NA NA NA 0.198 0.293 ...
    Number of loci: 2755
    Number of SNPs: 784 (28.46%)
    Number of heterozygous SNPs: 784 (100.00%)
    Chromosome: 1
    Segmenting DH signals...
     Segmenting by CBS...
      Chromosome: 1
       Random seed temporarily set (seed=c(10407, 1797822437, 388243314, 91406689, -519578635, -1381924756, 1089253019), kind="L'Ecuyer-CMRG")
     Segmenting by CBS...done
     List of 4
      $ data   :'data.frame':	2755 obs. of  4 variables:
       ..$ chromosome: int [1:2755] 1 1 1 1 1 1 1 1 1 1 ...
       ..$ x         : num [1:2755] 1.99e+08 1.99e+08 1.99e+08 1.99e+08 1.99e+08 ...
       ..$ y         : num [1:2755] NA NA NA 0.198 0.293 ...
       ..$ index     : int [1:2755] 1 2 3 4 5 6 7 8 9 10 ...
      $ output :'data.frame':	1 obs. of  6 variables:
       ..$ sampleName: chr NA
       ..$ chromosome: int 1
       ..$ start     : num 1.99e+08
       ..$ end       : num 2.07e+08
       ..$ nbrOfLoci : int 784
       ..$ mean      : num 0.22
      $ segRows:'data.frame':	1 obs. of  2 variables:
       ..$ startRow: int 4
       ..$ endRow  : int 2754
      $ params :List of 5
       ..$ alpha        : num 0.001
       ..$ undo         : num 0
       ..$ joinSegments : logi TRUE
       ..$ knownSegments:'data.frame':	1 obs. of  3 variables:
       .. ..$ chromosome: int 1
       .. ..$ start     : num 1.99e+08
       .. ..$ end       : num 2.07e+08
       ..$ seed         : int [1:7] 10407 1797822437 388243314 91406689 -519578635 -1381..
      - attr(*, "class")= chr [1:2] "CBS" "AbstractCBS"
      - attr(*, "processingTime")= 'proc_time' Named num [1:5] 0.009 0 0.009 0 0
       ..- attr(*, "names")= chr [1:5] "user.self" "sys.self" "elapsed" "user.child" ...
      - attr(*, "pkgDetails")= chr "DNAcopy v1.64.0"
      - attr(*, "randomSeed")= int [1:7] 10407 1797822437 388243314 91406689 -519578635 ..
     DH segmentation (locally-indexed) rows:
       startRow endRow
     1        4   2754
      int [1:2755] 54948 54949 54950 54951 54952 54953 54954 54955 54956 54957 ...
     DH segmentation rows:
       startRow endRow
     1    54951  57701
    Segmenting DH signals...done
    DH segmentation table:
       dhStart    dhEnd dhNbrOfLoci dhMean
    1 1.99e+08 2.07e+08         784   0.22
      startRow endRow
    1    54951  57701
    Rows:
    [1] 5
    TCN segmentation rows:
      startRow endRow
    5    54948  57702
    TCN and DH segmentation rows:
      startRow endRow
    5    54948  57702
      startRow endRow
    1    54951  57701
      startRow endRow
    1        1  37495
    2        1  37495
    3        1  37495
    4        1  37495
    5        1  37495
    6       NA     NA
    7    37496  50929
    8    50930  54947
    TCN segmentation (expanded) rows:
       startRow endRow
    1         1  37495
    2         1  37495
    3         1  37495
    4         1  37495
    5         1  37495
    6        NA     NA
    7     37496  50929
    8     50930  54947
    51    54948  57702
    TCN and DH segmentation rows:
      startRow endRow
    1        1  37495
    2       NA     NA
    3    37496  50929
    4    50930  54947
    5    54948  57702
    6    57703  57716
    7    57717  73297
      startRow endRow
    1        7    684
    2      693  12093
    3    12094  37250
    4    37253  37320
    5    37332  37449
    6       NA     NA
    7    37500  50929
    8    50930  54940
    9    54951  57701
      startRow endRow
    1        1  37495
    2        1  37495
    3        1  37495
    4        1  37495
    5        1  37495
    6       NA     NA
    7    37496  50929
    8    50930  54947
    9    54948  57702
    Total CN segmentation table (expanded):
      chromosome tcnStart   tcnEnd tcnNbrOfLoci tcnMean tcnNbrOfSNPs tcnNbrOfHets
    5          1 1.99e+08 2.07e+08         2755    2.59          784          784
    (TCN,DH) segmentation for one total CN segment:
      tcnId dhId chromosome tcnStart   tcnEnd tcnNbrOfLoci tcnMean tcnNbrOfSNPs
    5     5    1          1 1.99e+08 2.07e+08         2755    2.59          784
      tcnNbrOfHets  dhStart    dhEnd dhNbrOfLoci dhMean
    5          784 1.99e+08 2.07e+08         784   0.22
   Total CN segment #5 ([1.99122e+08,2.06513e+08]) of 7...done
   Total CN segment #6 ([2.06513e+08,2.06521e+08]) of 7...
    Number of TCN loci in segment: 14
    Locus data for TCN segment:
    'data.frame':	14 obs. of  9 variables:
     $ chromosome: int  1 1 1 1 1 1 1 1 1 1 ...
     $ x         : num  2.07e+08 2.07e+08 2.07e+08 2.07e+08 2.07e+08 ...
     $ CT        : num  3.55 4.92 3.73 3.45 3.18 ...
     $ betaT     : num  0.976 0.188 0.751 0.968 0.091 ...
     $ betaTN    : num  1.0103 0.2792 0.6263 1.015 0.0218 ...
     $ betaN     : num  0.9655 0.3359 0.6671 0.953 0.0692 ...
     $ muN       : num  1 0.5 0.5 1 0 0.5 0.5 0.5 0.5 0.5 ...
     $ index     : int  57703 57704 57705 57706 57707 57708 57709 57710 57711 57712 ...
     $ rho       : num  NA 0.442 0.253 NA NA ...
    Number of loci: 14
    Number of SNPs: 9 (64.29%)
    Number of heterozygous SNPs: 9 (100.00%)
    Chromosome: 1
    Segmenting DH signals...
     Segmenting by CBS...
      Chromosome: 1
       Random seed temporarily set (seed=c(10407, 1797822437, 388243314, 91406689, -519578635, -1381924756, 1089253019), kind="L'Ecuyer-CMRG")
     Segmenting by CBS...done
     List of 4
      $ data   :'data.frame':	14 obs. of  4 variables:
       ..$ chromosome: int [1:14] 1 1 1 1 1 1 1 1 1 1 ...
       ..$ x         : num [1:14] 2.07e+08 2.07e+08 2.07e+08 2.07e+08 2.07e+08 ...
       ..$ y         : num [1:14] NA 0.442 0.253 NA NA ...
       ..$ index     : int [1:14] 1 2 3 4 5 6 7 8 9 10 ...
      $ output :'data.frame':	1 obs. of  6 variables:
       ..$ sampleName: chr NA
       ..$ chromosome: int 1
       ..$ start     : num 2.07e+08
       ..$ end       : num 2.07e+08
       ..$ nbrOfLoci : int 9
       ..$ mean      : num 0.277
      $ segRows:'data.frame':	1 obs. of  2 variables:
       ..$ startRow: int 2
       ..$ endRow  : int 14
      $ params :List of 5
       ..$ alpha        : num 0.001
       ..$ undo         : num 0
       ..$ joinSegments : logi TRUE
       ..$ knownSegments:'data.frame':	1 obs. of  3 variables:
       .. ..$ chromosome: int 1
       .. ..$ start     : num 2.07e+08
       .. ..$ end       : num 2.07e+08
       ..$ seed         : int [1:7] 10407 1797822437 388243314 91406689 -519578635 -1381..
      - attr(*, "class")= chr [1:2] "CBS" "AbstractCBS"
      - attr(*, "processingTime")= 'proc_time' Named num [1:5] 0.002 0 0.002 0 0
       ..- attr(*, "names")= chr [1:5] "user.self" "sys.self" "elapsed" "user.child" ...
      - attr(*, "pkgDetails")= chr "DNAcopy v1.64.0"
      - attr(*, "randomSeed")= int [1:7] 10407 1797822437 388243314 91406689 -519578635 ..
     DH segmentation (locally-indexed) rows:
       startRow endRow
     1        2     14
      int [1:14] 57703 57704 57705 57706 57707 57708 57709 57710 57711 57712 ...
     DH segmentation rows:
       startRow endRow
     1    57704  57716
    Segmenting DH signals...done
    DH segmentation table:
       dhStart    dhEnd dhNbrOfLoci dhMean
    1 2.07e+08 2.07e+08           9  0.277
      startRow endRow
    1    57704  57716
    Rows:
    [1] 6
    TCN segmentation rows:
      startRow endRow
    6    57703  57716
    TCN and DH segmentation rows:
      startRow endRow
    6    57703  57716
      startRow endRow
    1    57704  57716
      startRow endRow
    1        1  37495
    2        1  37495
    3        1  37495
    4        1  37495
    5        1  37495
    6       NA     NA
    7    37496  50929
    8    50930  54947
    9    54948  57702
    TCN segmentation (expanded) rows:
       startRow endRow
    1         1  37495
    2         1  37495
    3         1  37495
    4         1  37495
    5         1  37495
    6        NA     NA
    7     37496  50929
    8     50930  54947
    9     54948  57702
    61    57703  57716
    TCN and DH segmentation rows:
      startRow endRow
    1        1  37495
    2       NA     NA
    3    37496  50929
    4    50930  54947
    5    54948  57702
    6    57703  57716
    7    57717  73297
       startRow endRow
    1         7    684
    2       693  12093
    3     12094  37250
    4     37253  37320
    5     37332  37449
    6        NA     NA
    7     37500  50929
    8     50930  54940
    9     54951  57701
    10    57704  57716
       startRow endRow
    1         1  37495
    2         1  37495
    3         1  37495
    4         1  37495
    5         1  37495
    6        NA     NA
    7     37496  50929
    8     50930  54947
    9     54948  57702
    10    57703  57716
    Total CN segmentation table (expanded):
      chromosome tcnStart   tcnEnd tcnNbrOfLoci tcnMean tcnNbrOfSNPs tcnNbrOfHets
    6          1 2.07e+08 2.07e+08           14    3.87            9            9
    (TCN,DH) segmentation for one total CN segment:
      tcnId dhId chromosome tcnStart   tcnEnd tcnNbrOfLoci tcnMean tcnNbrOfSNPs
    6     6    1          1 2.07e+08 2.07e+08           14    3.87            9
      tcnNbrOfHets  dhStart    dhEnd dhNbrOfLoci dhMean
    6            9 2.07e+08 2.07e+08           9  0.277
   Total CN segment #6 ([2.06513e+08,2.06521e+08]) of 7...done
   Total CN segment #7 ([2.06521e+08,2.47165e+08]) of 7...
    Number of TCN loci in segment: 15581
    Locus data for TCN segment:
    'data.frame':	15581 obs. of  9 variables:
     $ chromosome: int  1 1 1 1 1 1 1 1 1 1 ...
     $ x         : num  2.07e+08 2.07e+08 2.07e+08 2.07e+08 2.07e+08 ...
     $ CT        : num  3 3.06 2.61 2.76 2.21 ...
     $ betaT     : num  0.133 0.903 0.25 0.912 0.94 ...
     $ betaTN    : num  0.0458 1.0629 0.0301 1.1428 1.0005 ...
     $ betaN     : num  0.0875 0.8402 0.2201 0.7696 0.9395 ...
     $ muN       : num  0 1 0 1 1 0 1 1 1 0 ...
     $ index     : int  57717 57718 57719 57720 57721 57722 57723 57724 57725 57726 ...
     $ rho       : num  NA NA NA NA NA NA NA NA NA NA ...
    Number of loci: 15581
    Number of SNPs: 4492 (28.83%)
    Number of heterozygous SNPs: 4492 (100.00%)
    Chromosome: 1
    Segmenting DH signals...
     Segmenting by CBS...
      Chromosome: 1
       Random seed temporarily set (seed=c(10407, 1797822437, 388243314, 91406689, -519578635, -1381924756, 1089253019), kind="L'Ecuyer-CMRG")
     Segmenting by CBS...done
     List of 4
      $ data   :'data.frame':	15581 obs. of  4 variables:
       ..$ chromosome: int [1:15581] 1 1 1 1 1 1 1 1 1 1 ...
       ..$ x         : num [1:15581] 2.07e+08 2.07e+08 2.07e+08 2.07e+08 2.07e+08 ...
       ..$ y         : num [1:15581] NA NA NA NA NA NA NA NA NA NA ...
       ..$ index     : int [1:15581] 1 2 3 4 5 6 7 8 9 10 ...
      $ output :'data.frame':	1 obs. of  6 variables:
       ..$ sampleName: chr NA
       ..$ chromosome: int 1
       ..$ start     : num 2.07e+08
       ..$ end       : num 2.47e+08
       ..$ nbrOfLoci : int 4492
       ..$ mean      : num 0.229
      $ segRows:'data.frame':	1 obs. of  2 variables:
       ..$ startRow: int 30
       ..$ endRow  : int 15556
      $ params :List of 5
       ..$ alpha        : num 0.001
       ..$ undo         : num 0
       ..$ joinSegments : logi TRUE
       ..$ knownSegments:'data.frame':	1 obs. of  3 variables:
       .. ..$ chromosome: int 1
       .. ..$ start     : num 2.07e+08
       .. ..$ end       : num 2.47e+08
       ..$ seed         : int [1:7] 10407 1797822437 388243314 91406689 -519578635 -1381..
      - attr(*, "class")= chr [1:2] "CBS" "AbstractCBS"
      - attr(*, "processingTime")= 'proc_time' Named num [1:5] 0.035 0 0.035 0 0
       ..- attr(*, "names")= chr [1:5] "user.self" "sys.self" "elapsed" "user.child" ...
      - attr(*, "pkgDetails")= chr "DNAcopy v1.64.0"
      - attr(*, "randomSeed")= int [1:7] 10407 1797822437 388243314 91406689 -519578635 ..
     DH segmentation (locally-indexed) rows:
       startRow endRow
     1       30  15556
      int [1:15581] 57717 57718 57719 57720 57721 57722 57723 57724 57725 57726 ...
     DH segmentation rows:
       startRow endRow
     1    57746  73272
    Segmenting DH signals...done
    DH segmentation table:
       dhStart    dhEnd dhNbrOfLoci dhMean
    1 2.07e+08 2.47e+08        4492  0.229
      startRow endRow
    1    57746  73272
    Rows:
    [1] 7
    TCN segmentation rows:
      startRow endRow
    7    57717  73297
    TCN and DH segmentation rows:
      startRow endRow
    7    57717  73297
      startRow endRow
    1    57746  73272
       startRow endRow
    1         1  37495
    2         1  37495
    3         1  37495
    4         1  37495
    5         1  37495
    6        NA     NA
    7     37496  50929
    8     50930  54947
    9     54948  57702
    10    57703  57716
    TCN segmentation (expanded) rows:
       startRow endRow
    1         1  37495
    2         1  37495
    3         1  37495
    4         1  37495
    5         1  37495
    6        NA     NA
    7     37496  50929
    8     50930  54947
    9     54948  57702
    10    57703  57716
    71    57717  73297
    TCN and DH segmentation rows:
      startRow endRow
    1        1  37495
    2       NA     NA
    3    37496  50929
    4    50930  54947
    5    54948  57702
    6    57703  57716
    7    57717  73297
       startRow endRow
    1         7    684
    2       693  12093
    3     12094  37250
    4     37253  37320
    5     37332  37449
    6        NA     NA
    7     37500  50929
    8     50930  54940
    9     54951  57701
    10    57704  57716
    11    57746  73272
       startRow endRow
    1         1  37495
    2         1  37495
    3         1  37495
    4         1  37495
    5         1  37495
    6        NA     NA
    7     37496  50929
    8     50930  54947
    9     54948  57702
    10    57703  57716
    11    57717  73297
    Total CN segmentation table (expanded):
      chromosome tcnStart   tcnEnd tcnNbrOfLoci tcnMean tcnNbrOfSNPs tcnNbrOfHets
    7          1 2.07e+08 2.47e+08        15581    2.64         4492         4492
    (TCN,DH) segmentation for one total CN segment:
      tcnId dhId chromosome tcnStart   tcnEnd tcnNbrOfLoci tcnMean tcnNbrOfSNPs
    7     7    1          1 2.07e+08 2.47e+08        15581    2.64         4492
      tcnNbrOfHets  dhStart    dhEnd dhNbrOfLoci dhMean
    7         4492 2.07e+08 2.47e+08        4492  0.229
   Total CN segment #7 ([2.06521e+08,2.47165e+08]) of 7...done
      chromosome tcnId dhId tcnStart   tcnEnd tcnNbrOfLoci tcnMean tcnNbrOfSNPs
   1           1     1    1 5.54e+05 1.21e+08        37495    1.38        10146
   2           1     1    2 5.54e+05 1.21e+08        37495    1.38        10146
   3           1     1    3 5.54e+05 1.21e+08        37495    1.38        10146
   4           1     1    4 5.54e+05 1.21e+08        37495    1.38        10146
   5           1     1    5 5.54e+05 1.21e+08        37495    1.38        10146
   6           1     2    1 1.21e+08 1.42e+08            0      NA            0
   7           1     3    1 1.42e+08 1.86e+08        13434    2.07         3770
   8           1     4    1 1.86e+08 1.99e+08         4018    2.71         1271
   9           1     5    1 1.99e+08 2.07e+08         2755    2.59          784
   10          1     6    1 2.07e+08 2.07e+08           14    3.87            9
   11          1     7    1 2.07e+08 2.47e+08        15581    2.64         4492
      tcnNbrOfHets  dhStart    dhEnd dhNbrOfLoci dhMean
   1         10146 5.54e+05 4.28e+06         169 0.4145
   2         10146 4.28e+06 4.27e+07        3260 0.4944
   3         10146 4.27e+07 1.20e+08        6657 0.5205
   4         10146 1.20e+08 1.20e+08           8 0.0767
   5         10146 1.20e+08 1.21e+08          52 0.5123
   6             0       NA       NA          NA     NA
   7          3770 1.42e+08 1.86e+08        3770 0.0943
   8          1271 1.86e+08 1.99e+08        1271 0.2563
   9           784 1.99e+08 2.07e+08         784 0.2197
   10            9 2.07e+08 2.07e+08           9 0.2769
   11         4492 2.07e+08 2.47e+08        4492 0.2290
   Calculating (C1,C2) per segment...
   Calculating (C1,C2) per segment...done
   Number of segments: 11
  Segmenting paired tumor-normal signals using Paired PSCBS...done
  Post-segmenting TCNs...
   Number of segments: 11
   Number of chromosomes: 1
   [1] 1
   Chromosome 1 ('chr01') of 1...
    Rows:
     [1]  1  2  3  4  5  6  7  8  9 10 11
    Number of segments: 11
    TCN segment #1 ('1') of 7...
     Rows:
     [1] 1 2 3 4 5
     TCN & DH segRows before:
       tcn.startRow tcn.endRow dh.startRow dh.endRow
     1            1      37495           7       684
     2            1      37495         693     12093
     3            1      37495       12094     37250
     4            1      37495       37253     37320
     5            1      37495       37332     37449
     Range [1,37495]
     DH segment #1 of 5...
      [xStart,xEnd] = [554484,4278527]
      [idxStart,idxEnd] = [1,687]
      Number of TCN loci: 687
      [idxStart,idxEnd] = [1,687]
      Number of non-missing TCN loci: 687
     DH segment #1 of 5...done
     DH segment #2 of 5...
      [xStart,xEnd] = [4278527,42651414]
      [idxStart,idxEnd] = [688,12093]
      Number of TCN loci: 11406
      [idxStart,idxEnd] = [688,12093]
      Number of non-missing TCN loci: 11406
     DH segment #2 of 5...done
     DH segment #3 of 5...
      [xStart,xEnd] = [42651414,119796080]
      [idxStart,idxEnd] = [12094,37252]
      Number of TCN loci: 25159
      [idxStart,idxEnd] = [12094,37252]
      Number of non-missing TCN loci: 25159
     DH segment #3 of 5...done
     DH segment #4 of 5...
      [xStart,xEnd] = [119796080,119932126]
      [idxStart,idxEnd] = [37253,37324]
      Number of TCN loci: 72
      [idxStart,idxEnd] = [37253,37324]
      Number of non-missing TCN loci: 72
     DH segment #4 of 5...done
     DH segment #5 of 5...
      [xStart,xEnd] = [119932126,120992603]
      [idxStart,idxEnd] = [37325,37495]
      Number of TCN loci: 171
      [idxStart,idxEnd] = [37325,37495]
      Number of non-missing TCN loci: 171
     DH segment #5 of 5...done
     TCN & DH segRows afterward:
       tcn.startRow tcn.endRow dh.startRow dh.endRow
     1            1        687           7       684
     2          688      12093         693     12093
     3        12094      37252       12094     37250
     4        37253      37324       37253     37320
     5        37325      37495       37332     37449
     Number of TCNs before: 37495
     Number of TCNs after: 37495
    TCN segment #1 ('1') of 7...done
    TCN segment #2 ('2') of 7...
     Nothing todo. Only one DH segmentation. Skipping.
    TCN segment #2 ('2') of 7...done
    TCN segment #3 ('3') of 7...
     Nothing todo. Only one DH segmentation. Skipping.
    TCN segment #3 ('3') of 7...done
    TCN segment #4 ('4') of 7...
     Nothing todo. Only one DH segmentation. Skipping.
    TCN segment #4 ('4') of 7...done
    TCN segment #5 ('5') of 7...
     Nothing todo. Only one DH segmentation. Skipping.
    TCN segment #5 ('5') of 7...done
    TCN segment #6 ('6') of 7...
     Nothing todo. Only one DH segmentation. Skipping.
    TCN segment #6 ('6') of 7...done
    TCN segment #7 ('7') of 7...
     Nothing todo. Only one DH segmentation. Skipping.
    TCN segment #7 ('7') of 7...done
   Chromosome 1 ('chr01') of 1...done
   Update (C1,C2) per segment...
   Update (C1,C2) per segment...done
  Post-segmenting TCNs...done
    chromosome tcnId dhId tcnStart   tcnEnd tcnNbrOfLoci tcnMean tcnNbrOfSNPs
  1          1     1    1 5.54e+05 4.28e+06          687    1.40          169
  2          1     1    2 4.28e+06 4.27e+07        11406    1.38         3260
  3          1     1    3 4.27e+07 1.20e+08        25159    1.38         6657
  4          1     1    4 1.20e+08 1.20e+08           72    1.47            8
  5          1     1    5 1.20e+08 1.21e+08          171    1.44           52
  6          1     2    1 1.21e+08 1.42e+08            0      NA            0
    tcnNbrOfHets  dhStart    dhEnd dhNbrOfLoci dhMean c1Mean c2Mean
  1          169 5.54e+05 4.28e+06         169 0.4145  0.411  0.992
  2         3260 4.28e+06 4.27e+07        3260 0.4944  0.348  1.030
  3         6657 4.27e+07 1.20e+08        6657 0.5205  0.332  1.052
  4            8 1.20e+08 1.20e+08           8 0.0767  0.679  0.792
  5           52 1.20e+08 1.21e+08          52 0.5123  0.351  1.089
  6            0       NA       NA          NA     NA     NA     NA
     chromosome tcnId dhId tcnStart   tcnEnd tcnNbrOfLoci tcnMean tcnNbrOfSNPs
  6           1     2    1 1.21e+08 1.42e+08            0      NA            0
  7           1     3    1 1.42e+08 1.86e+08        13434    2.07         3770
  8           1     4    1 1.86e+08 1.99e+08         4018    2.71         1271
  9           1     5    1 1.99e+08 2.07e+08         2755    2.59          784
  10          1     6    1 2.07e+08 2.07e+08           14    3.87            9
  11          1     7    1 2.07e+08 2.47e+08        15581    2.64         4492
     tcnNbrOfHets  dhStart    dhEnd dhNbrOfLoci dhMean c1Mean c2Mean
  6             0       NA       NA          NA     NA     NA     NA
  7          3770 1.42e+08 1.86e+08        3770 0.0943  0.935   1.13
  8          1271 1.86e+08 1.99e+08        1271 0.2563  1.007   1.70
  9           784 1.99e+08 2.07e+08         784 0.2197  1.009   1.58
  10            9 2.07e+08 2.07e+08           9 0.2769  1.400   2.47
  11         4492 2.07e+08 2.47e+08        4492 0.2290  1.017   1.62
  Calling ROH...
   Segment #1 of 11...
    Calling ROH for a single segment...
     Number of SNPs: 687
    Calling ROH for a single segment...done
   Segment #1 of 11...done
   Segment #2 of 11...
    Calling ROH for a single segment...
     Number of SNPs: 11406
    Calling ROH for a single segment...done
   Segment #2 of 11...done
   Segment #3 of 11...
    Calling ROH for a single segment...
     Number of SNPs: 25159
    Calling ROH for a single segment...done
   Segment #3 of 11...done
   Segment #4 of 11...
    Calling ROH for a single segment...
     Number of SNPs: 72
    Calling ROH for a single segment...done
   Segment #4 of 11...done
   Segment #5 of 11...
    Calling ROH for a single segment...
     Number of SNPs: 171
    Calling ROH for a single segment...done
   Segment #5 of 11...done
   Segment #6 of 11...
    Calling ROH for a single segment...
     Number of SNPs: 0
    Calling ROH for a single segment...done
   Segment #6 of 11...done
   Segment #7 of 11...
    Calling ROH for a single segment...
     Number of SNPs: 13434
    Calling ROH for a single segment...done
   Segment #7 of 11...done
   Segment #8 of 11...
    Calling ROH for a single segment...
     Number of SNPs: 4018
    Calling ROH for a single segment...done
   Segment #8 of 11...done
   Segment #9 of 11...
    Calling ROH for a single segment...
     Number of SNPs: 2755
    Calling ROH for a single segment...done
   Segment #9 of 11...done
   Segment #10 of 11...
    Calling ROH for a single segment...
     Number of SNPs: 14
    Calling ROH for a single segment...done
   Segment #10 of 11...done
   Segment #11 of 11...
    Calling ROH for a single segment...
     Number of SNPs: 15581
    Calling ROH for a single segment...done
   Segment #11 of 11...done
   ROH calls:
    logi [1:11] FALSE FALSE FALSE TRUE FALSE NA ...
      Mode   FALSE    TRUE    NA's 
   logical       9       1       1 
  Calling ROH...done
  Error: processing vignette 'PairedPSCBS.tex.rsp' failed with diagnostics:
  there is no package called ‘Hmisc’
  --- failed re-building ‘PairedPSCBS.tex.rsp’
  
  SUMMARY: processing the following file failed:
    ‘PairedPSCBS.tex.rsp’
  
  Error: Vignette re-building failed.
  Execution halted

2 errors ✖ | 1 warning ✖ | 3 notes ✖

── PSCBS 0.66.0: OK

  Build ID:   PSCBS_0.66.0.tar.gz-f32c847c506f4ccca0736e356b6eb7bf
  Platform:   macOS 10.13.6 High Sierra, R-release, CRAN's setup
  Submitted:  27m 57.2s ago
  Build time: 13m 22.9s

0 errors ✔ | 0 warnings ✔ | 0 notes ✔

── PSCBS 0.66.0: ERROR

  Build ID:   PSCBS_0.66.0.tar.gz-5103c1b9bc9d4f1fa4053c8baddfd025
  Platform:   Apple Silicon (M1), macOS 11.6 Big Sur, R-release
  Submitted:  27m 57.2s ago
  Build time: 3m 29.1s

❯ checking package dependencies ... ERROR
  Package suggested but not available: ‘Hmisc’
  
  The suggested packages are required for a complete check.
  Checking can be attempted without them by setting the environment
  variable _R_CHECK_FORCE_SUGGESTS_ to a false value.
  
  See section ‘The DESCRIPTION file’ in the ‘Writing R Extensions’
  manual.

1 error ✖ | 0 warnings ✔ | 0 notes ✔

── PSCBS 0.66.0: OK

  Build ID:   PSCBS_0.66.0.tar.gz-324cf49ae8414efa88b5cc3455484f78
  Platform:   Oracle Solaris 10, x86, 32 bit, R release, Oracle Developer Studio 12.6
  Submitted:  27m 57.2s ago
  Build time: 20m 23.1s

0 errors ✔ | 0 warnings ✔ | 0 notes ✔

── PSCBS 0.66.0: ERROR

  Build ID:   PSCBS_0.66.0.tar.gz-a82d7b7c445340fc8dc998517d363db6
  Platform:   Windows Server 2008 R2 SP1, R-devel, 32/64 bit
  Submitted:  27m 57.2s ago
  Build time: 1m 58.8s

❯ checking package dependencies ... ERROR
  Package required but not available: 'DNAcopy'
  
  See section 'The DESCRIPTION file' in the 'Writing R Extensions'
  manual.

1 error ✖ | 0 warnings ✔ | 0 notes ✔

── PSCBS 0.66.0: OK

  Build ID:   PSCBS_0.66.0.tar.gz-c26455fe7c1244d0824fbb7f59d796ac
  Platform:   Windows Server 2008 R2 SP1, R-release, 32/64 bit
  Submitted:  27m 57.2s ago
  Build time: 9m 5.5s

0 errors ✔ | 0 warnings ✔ | 0 notes ✔
> 
```
