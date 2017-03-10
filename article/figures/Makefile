SVGFILES = $(wildcard *.svg)
PDFFILES = $(SVGFILES:.svg=.pdf)
 
all: $(PDFFILES)

%.pdf : %.svg
	@echo "======= Converting SVG to PDF: $< ======="
	@inkscape --without-gui --file=$< --export-pdf=$@ --export-text-to-path
