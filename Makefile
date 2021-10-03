BIB = bibtex
PDF = pdflatex

all : main.pdf

FIGURES = $(wildcard ./PhD-text/*/*.jpg ./PhD-text/*/*.jpeg ./PhD-text/*/*.pdf ./PhD-text/*/*.png ./PhD-text/*/*.eps)

# It is important to run pdflatex twice to update thre refferences.
main.pdf : main.tex bibl.bib ./PhD-text/*/*.tex $(FIGURES)
	$(PDF) main.tex
	$(BIB) main.aux
	$(PDF) main.tex
	$(PDF) main.tex
	
clean:
	rm -f *.bbl
	rm -f *.blg
	rm -f *.log
	rm -f *.aux
	rm -f main.pdf
	rm -f *.dvi
	rm -f *.toc
	rm -f *.aux
	rm -f *.out

	
