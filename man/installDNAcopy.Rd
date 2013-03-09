%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  installDNAcopy.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{installDNAcopy}
\alias{installDNAcopy.default}
\alias{installDNAcopy}

\title{Install the DNAcopy package}

\usage{
\method{installDNAcopy}{default}(..., force=FALSE)
}

\description{
  Install the DNAcopy package, if missing.
}

\arguments{
  \item{...}{Arguments passed to the install function.}
  \item{force}{If \code{\link[base:logical]{FALSE}} and the \pkg{DNAcopy} package is already
    installed, then it will not be re-install.
    If \code{\link[base:logical]{TRUE}}, it will be installed.}
}

\value{
  Returns nothing.
}

\details{
  This function is will download and call the \code{biocLite()}
  installation function from the Bioconductor Project website.
  This function will also make sure that \pkg{DNAcopy} is loaded so
  that it is reported by \code{\link[utils]{sessionInfo}}.
}

\author{Henrik Bengtsson}


\keyword{internal}
