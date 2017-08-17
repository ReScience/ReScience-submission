all: stojic2017.tex stojic2017.pdf clean

stojic2017.tex: stojic2017.md	
	pandoc --standalone --filter ~/.cabal/bin/pandoc-crossref --template=rescience-template.tex --latex-engine=xelatex --biblatex --bibliography=stojic2017.bib -M "crossrefYaml=crossref.yaml" --output stojic2017.tex stojic2017.md

stojic2017.pdf: stojic2017.tex
	xelatex stojic2017
	biber stojic2017
	xelatex stojic2017
	xelatex stojic2017

clean:
	rm -f *.aux *.out *.log *.nav *.toc *.snm *.ttt *.blg *.bbl *.fff *.suppinfo *latexmk *fls *synctex.gz *.bcf *.xml
