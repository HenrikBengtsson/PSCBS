<% if (git_branch() == "develop") { %>
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
future vignette '[A Future for R: Controlling Default Future Strategy](https://cran.r-project.org/web/packages/future/vignettes/future-5-startup.html)'.
<% } %>


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

<% } # if (git_branch() != "master") %>
