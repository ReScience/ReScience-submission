A reference implementation of
"*Fast-Activating Voltage- and Calcium-Dependent Potassium (BK) Conductance
Promotes Bursting in Pituitary Cells: A Dynamic Clamp Study*,
J. Tabak, M. Tomaiuolo, A. Gonzalez-Iglesias,  L. Milescu and R. Bertram,
Journal of Neuroscience 31.46 (2011), 10.1523/JNEUROSCI.3235-11.2011"

### Required tools for producing the pdf

You'll need [pandoc](http://pandoc.org) (a universal document converter) and a
full [TeX distribution](https://www.tug.org/texlive/).

For pandoc, you'll also need the
[pandoc-crossref](https://github.com/lierdakil/pandoc-crossref) filter that you can
easily install with:

```
$ cabal update
$ cabal install pandoc-crossref
```

### How to build the PDF?

In a console, type:

```
pandoc --standalone --filter ~/.cabal/bin/pandoc-crossref --template=rescience-template.tex --latex-engine=xelatex --biblatex --bibliography=bibliography.bib -M "crossrefYaml=crossref.yaml" --output tennøe-hodne-haug-weltzien-einevoll-halnes-2018.tex tennøe-hodne-haug-weltzien-einevoll-halnes-2018.md
xelatex tennøe-hodne-haug-weltzien-einevoll-halnes-2018
biber tennøe-hodne-haug-weltzien-einevoll-halnes-2018
xelatex tennøe-hodne-haug-weltzien-einevoll-halnes-2018
xelatex tennøe-hodne-haug-weltzien-einevoll-halnes-2018
```

Alternatively, you can also type `make`.
