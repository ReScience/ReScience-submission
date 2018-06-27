#!/bin/sh
pandoc \
  --standalone \
  --filter /usr/bin/pandoc-crossref \
  --latex-engine=xelatex \
  --template=rescience-template.tex \
  --biblatex --bibliography=bibliography.bib \
  -M "crossrefYaml=/crossref.yaml" \
  --output Hathway-Goodman-2018.tex Hathway-Goodman-2018.md
xelatex Hathway-Goodman-2018.tex >make.log
biber Hathway-Goodman-2018 >>make.log
xelatex Hathway-Goodman-2018.tex >>make.log
xelatex Hathway-Goodman-2018.tex >>make.log
