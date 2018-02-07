pandoc --standalone --filter pandoc-fignos --filter pandoc-tablenos --filter pandoc-eqnos --template=rescience-template.tex --latex-engine=xelatex --biblatex --bibliography=senden-schuecker-hahne-diesmann-goebel-2018.bib -M "crossrefYaml=crossref.yaml" --output senden-schuecker-hahne-diesmann-goebel-2018.tex senden-schuecker-hahne-diesmann-goebel-2018.md

xelatex senden-schuecker-hahne-diesmann-goebel-2018
biber senden-schuecker-hahne-diesmann-goebel-2018
xelatex senden-schuecker-hahne-diesmann-goebel-2018
xelatex senden-schuecker-hahne-diesmann-goebel-2018
