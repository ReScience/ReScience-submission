#!/bin/sh
docker run \
  -v `pwd`/:/root/ \
  pandoc \
  --standalone --latex-engine=xelatex \
  --template=rescience-template.tex \
  --biblatex --bibliography=bibliography.bib \
  -M "crossrefYaml=/crossref.yaml" \
  --output Hathway-Goodman-2018.pdf Hathway-Goodman-2018.md
