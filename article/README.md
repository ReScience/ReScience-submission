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

```
cd ~/work/0\ research/git\ projects/ReScience-submission/article
pandoc --standalone --filter /Users/owenpetchey/Library/Haskell/bin/pandoc-crossref --template=rescience-template.tex --latex-engine=xelatex --biblatex --bibliography=article.bib -M "crossrefYaml=crossref.yaml" --output article.tex article.md
xelatex article
biber article
xelatex article
xelatex article

```

