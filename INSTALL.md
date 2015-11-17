## Installation
R package <%=pkg()%> is available on [CRAN](http://cran.r-project.org/package=<%=pkg()%>) and can be installed in R as:
```r
install.packages('<%=pkg()%>')
```

<% if (git_branch() != "master" && !grepl("^release/", git_branch())) { %>
### Pre-release version

To install the pre-release version that is available in branch `<%=git_branch()%>`, use:
```r
source('http://callr.org/install#<%=github_repos()%>@<%=git_branch()%>')
```
This will install the package from source.  <% if (file.exists("src")) { %><%-%>
Because of this and because this package also compiles native code,
Windows users need to have
[Rtools](https://cran.r-project.org/bin/windows/Rtools/) installed and
OS X users need to have [Xcode](https://developer.apple.com/xcode/)
installed.
<% } # if (file.exists("src")) %>

<% if (git_branch() == "feature/future") { %>
#### Parallel processing
The `<%=git_branch()%>` branch supports segmentation of the
chromosomes in parallel (asynchronously) by adding the following
```r
future::plan("multicore")
```
to the beginning of the PSCBS script.  Everything else will work the
same.  To reset to non-parallel (synchronously) processing, use
`future::plan("eager")`.

An alternative to editing the R script is to set environment variable
`R_FUTURE_PLAN`, e.g.
```sh
export R_FUTURE_PLAN=multicore
```
To control the maximum number of cores the multicore processing may
use set environment variable `MC_CORES`, e.g.
```sh
export MC_CORES=4
```
This variable is defined and read by the 'parallel' package when it
is loaded (so not when R itself is started) and used to set options
`mc.cores`, which is acknowledged by `future::plan("multicore")`.

This above also be set in the cross-platform `~/.Renviron` file as:
```r
R_FUTURE_PLAN=multicore
MC_CORES=4
```

<% } %>

<% } # if (git_branch() != "master") %>
