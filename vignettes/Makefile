# Makefile to use knitr for package vignettes
# (modified from knitr package source)
#
# DO NOT EXCLUDE FROM BUILD!

PDFS = ResponseModels-knitr.pdf

all: $(PDFS)

clean:
	rm -rf *-concordance.tex *.bbl *.blg *.aux *.out *.log *.spl *.md figure/ *.idx *.ilg *.ind *.toc

#; tools::compactPDF('$*.pdf', gs_quality = 'ebook')
%.pdf: %.Rnw
	$(R_HOME)/bin/Rscript -e "if (getRversion() < '3.0.0') knitr::knit2pdf('$*.Rnw') else tools::texi2pdf('$*.tex')"

%.html: %.Rmd
	$(R_HOME)/bin/Rscript -e "if (getRversion() < '3.0.0') knitr::knit2html('$*.Rmd')"
