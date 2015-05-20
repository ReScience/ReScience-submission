### ReScience submission repository

This is the submission repository for the [Re**Science** journal](https://github.com/ReScience/ReScience/wiki).

### How to build the PDF ?

In a console, type:

```
pandoc --standalone --template=ReScience-template.tex --latex-engine=xelatex --biblatex --bibliography=article.bib --output article.tex article.md
xelatex article
biber article
xelatex article
xelatex article
```
