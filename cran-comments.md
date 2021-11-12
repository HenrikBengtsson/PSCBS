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
| devel         | L M    | L          |    W            |

*Legend: OS: L = Linux, S = Solaris, M = macOS, M1 = macOS M1, W = Windows*


R-hub checks:

```
> res <- rhub::check(platform = c(
  "debian-clang-devel", "debian-gcc-patched", "linux-x86_64-centos-epel",
  "macos-highsierra-release-cran", "solaris-x86-patched-ods",
  "windows-x86_64-release"))
> res

── PSCBS 0.65.0-9002: OK

  Build ID:   PSCBS_0.65.0-9002.tar.gz-355deceb5d1b45239b61795ae37279f1
  Platform:   Debian Linux, R-devel, clang, ISO-8859-15 locale
  Submitted:  25m 17.9s ago
  Build time: 16m 20.2s

0 errors ✔ | 0 warnings ✔ | 0 notes ✔

── PSCBS 0.65.0-9002: OK

  Build ID:   PSCBS_0.65.0-9002.tar.gz-74092b87421e4660b97884b9da623140
  Platform:   Debian Linux, R-patched, GCC
  Submitted:  25m 17.9s ago
  Build time: 14m 2.5s

0 errors ✔ | 0 warnings ✔ | 0 notes ✔

── PSCBS 0.65.0-9002: NOTE

  Build ID:   PSCBS_0.65.0-9002.tar.gz-6590c41d4e7a40d191c33e9edbb1f150
  Platform:   CentOS 8, stock R from EPEL
  Submitted:  25m 18s ago
  Build time: 10m 48.4s

❯ checking package dependencies ... NOTE
  Package suggested but not available for checking: ‘Hmisc’

0 errors ✔ | 0 warnings ✔ | 1 note ✖

── PSCBS 0.65.0-9002: OK

  Build ID:   PSCBS_0.65.0-9002.tar.gz-7e50c30cb4d74188a99e2dd9653b4193
  Platform:   macOS 10.13.6 High Sierra, R-release, CRAN's setup
  Submitted:  25m 18s ago
  Build time: 13m 10.9s

0 errors ✔ | 0 warnings ✔ | 0 notes ✔

── PSCBS 0.65.0-9002: OK

  Build ID:   PSCBS_0.65.0-9002.tar.gz-c2778bb9a43f4516bebb8704df1445e8
  Platform:   Oracle Solaris 10, x86, 32 bit, R release, Oracle Developer Studio 12.6
  Submitted:  25m 18s ago
  Build time: 20m 9.7s

0 errors ✔ | 0 warnings ✔ | 0 notes ✔

── PSCBS 0.65.0-9002: OK

  Build ID:   PSCBS_0.65.0-9002.tar.gz-00112bd7405244eebaa1260bb3df840e
  Platform:   Windows Server 2008 R2 SP1, R-release, 32/64 bit
  Submitted:  25m 18s ago
  Build time: 8m 32.2s

0 errors ✔ | 0 warnings ✔ | 0 notes ✔
```

Special care:

```r
> res <- rhub::check(platform = "macos-m1-bigsur-release", env_vars = c("_R_CHECK_FORCE_SUGGESTS_" = "false"))
> res

── PSCBS 0.66.0: NOTE

  Build ID:   PSCBS_0.66.0.tar.gz-46d044ada3d74f1eba3a184754d0daf3
  Platform:   Apple Silicon (M1), macOS 11.6 Big Sur, R-release
  Submitted:  5m 31.8s ago
  Build time: 5m 27.9s

❯ checking package dependencies ... NOTE
  Package suggested but not available for checking: ‘Hmisc’

0 errors ✔ | 0 warnings ✔ | 1 note ✖
```
