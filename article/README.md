### How to build the PDF ?

In a console, type:

```
pandoc --standalone --template=ReScience-template.tex --latex-engine=xelatex --biblatex --bibliography=article.bib --output article.tex article.md
xelatex article
biber article
xelatex article
xelatex article
```
