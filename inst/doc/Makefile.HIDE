# MACRO DEFINITIONS
RM=rm -f
MV=mv -f
CP=cp -p

# TARGET DEFINITIONS
all: rnw index.html

rnw: *.Rnw
	R --vanilla -e "R.rsp::buildNonSweaveVignettes()"

index.html:
	R --vanilla -e "R.rsp::buildPkgIndexHtml()"

clean:
	$(RM) index.html.rsp
	$(RM) *.aux *.bbl *.blg *.log *.toc *.tex *.rsp.R
	$(RM) -r figures/
#	$(RM) Makefile
