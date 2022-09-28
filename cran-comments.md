# CRAN submission PSCBS 0.67.0

on 2022-09-27

Thanks in advance


## Notes not sent to CRAN

### R CMD check validation

The package has been verified using `R CMD check --as-cran` on:

| R version     | GitHub | R-hub    | mac/win-builder |
| ------------- | ------ | -------- | --------------- |
| 3.6.x         | L      |          |                 |
| 4.0.x         | L      |          |                 |
| 4.1.x         | L      |          |                 |
| 4.2.x         | L M W  | L M M1 W | M1 w            |
| devel         | L M W  | L        |    w            |

_Legend: OS: L = Linux, M = macOS, M1 = macOS M1, W = Windows_


R-hub checks:

```r
res <- rhub::check(platforms = c(
  "debian-clang-devel", "debian-gcc-patched",
  "macos-highsierra-release-cran", "macos-m1-bigsur-release",
  "windows-x86_64-release"))
print(res)
```

gives

```

── PSCBS 0.66.0-9004: OK

  Build ID:   PSCBS_0.66.0-9004.tar.gz-cc87215cc60040fcaadaf0c0bf496791
  Platform:   Debian Linux, R-devel, clang, ISO-8859-15 locale
  Submitted:  22m 20.8s ago
  Build time: 21m 38.2s

0 errors ✔ | 0 warnings ✔ | 0 notes ✔

── PSCBS 0.66.0-9004: PREPERROR

  Build ID:   PSCBS_0.66.0-9004.tar.gz-c2202a4ffbee460d805ee1da16a107f9
  Platform:   Debian Linux, R-patched, GCC
  Submitted:  22m 20.8s ago
  Build time: 20m 11.7s

❯ Build failed during preparation or aborted
...

── PSCBS 0.66.0-9004: OK

  Build ID:   PSCBS_0.66.0-9004.tar.gz-2940e83df25243858e2313512d3110ba
  Platform:   macOS 10.13.6 High Sierra, R-release, CRAN's setup
  Submitted:  22m 20.8s ago
  Build time: 5m 44.4s

0 errors ✔ | 0 warnings ✔ | 0 notes ✔

── PSCBS 0.66.0-9004: CREATED

  Build ID:   PSCBS_0.66.0-9004.tar.gz-b84a6bb86c124617a29257837353b4e4
  Platform:   Apple Silicon (M1), macOS 11.6 Big Sur, R-release
  Submitted:  22m 20.9s ago


── PSCBS 0.66.0-9004: OK

  Build ID:   PSCBS_0.66.0-9004.tar.gz-ebfbab7bf6ba43d1a18bb8cb57ae1dde
  Platform:   Windows Server 2022, R-release, 32/64 bit
  Submitted:  22m 20.9s ago
  Build time: 10m 21.4s

0 errors ✔ | 0 warnings ✔ | 0 notes ✔
```
