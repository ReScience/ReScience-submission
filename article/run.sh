pandoc --standalone --filter ~/.cabal/bin/pandoc-crossref --template=rescience-template.tex --latex-engine=xelatex --biblatex --bibliography=article.bib -M "crossrefYaml=crossref.yaml" --output article.tex article.md
xelatex article
biber article
xelatex article
xelatex article
