LATEX=pdflatex
LATEXOPT=--shell-escape

LATEXMK=latexmk
LATEXMKOPT=-pdf

MAIN=hf
HEADER=header
MAINDIR:=$(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))

all:
	mkdir -p compile
	cp $(MAIN).tex compile
	cp $(HEADER).tex compile
	cd compile && \
	$(LATEXMK) $(LATEXMKOPT) $(MAIN).tex
	cp $(MAINDIR)/compile/$(MAIN).pdf $(MAINDIR)
clean:
	rm -rf compile
veryclean:
	rm -rf compile
	rm $(MAIN).pdf

.PHONY: clean veryclean all
