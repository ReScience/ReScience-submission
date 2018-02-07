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

### How to build the PDF ?

In a console, type:

```shell
$ pandoc --standalone --filter pandoc-fignos --filter pandoc-tablenos --filter pandoc-eqnos --template=rescience-template.tex --latex-engine=xelatex --biblatex --bibliography=senden-schuecker-hahne-diesmann-goebel-2018.bib -M "crossrefYaml=crossref.yaml" --output senden-schuecker-hahne-diesmann-goebel-2018.tex senden-schuecker-hahne-diesmann-goebel-2018.md

$ xelatex senden-schuecker-hahne-diesmann-goebel-2018
$ biber senden-schuecker-hahne-diesmann-goebel-2018
$ xelatex senden-schuecker-hahne-diesmann-goebel-2018
$ xelatex senden-schuecker-hahne-diesmann-goebel-2018
```

Alternativaley, you can also run `compile_pdf.sh`.
