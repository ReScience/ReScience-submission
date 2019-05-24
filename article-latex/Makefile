# From  https://tex.stackexchange.com/questions/40738/how-to-properly-make-a-latex-project
# -----------------------------------------------------------------------------
# You want latexmk to *always* run, because make does not have all the info.
# Also, include non-file targets in .PHONY so they are run regardless of any
# file of the given name existing.
.PHONY: clean

# The first rule in a Makefile is the one executed by default ("make"). It
# should always be the "all" rule, so that "make" and "make all" are identical.

all: article
article: article.pdf


# CUSTOM BUILD RULES
# -----------------------------------------------------------------------------
metadata.tex: metadata.yaml
	./yaml-to-latex.py -i $< -o $@


# MAIN LATEXMK RULE
# -----------------------------------------------------------------------------
# -pdf tells latexmk to generate PDF directly (instead of DVI).
# -pdflatex="" tells latexmk to call a specific backend with specific options.
# -use-make tells latexmk to call make for generating missing files.
# -interaction=nonstopmode keeps the pdflatex backend from stopping at a
# missing file reference and interactively asking you for an alternative.

article.pdf: article.tex content.tex metadata.tex rescience.cls
	latexmk -pdf -pdflatex="xelatex -interaction=nonstopmode" -use-make article.tex

clean:
	@latexmk -CA
	@rm -f *.bbl
	@rm -f *.run.xml
