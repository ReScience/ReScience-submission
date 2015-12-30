images:
	convert print-mesh-1.png -gravity South -resize x340\! print-mesh-2.png
	convert print-mesh-2.png -gravity South -extent 472x398 print-mesh-3.png
	montage coord-id-1.png print-mesh-3.png -tile 2x -geometry +1+1 -bordercolor white -border 5x5 node-edge_selection.png
	montage acc_cost-1.png impute_na-1.png -gravity North -tile 2x -geometry +1+1 -bordercolor white -border 5x5 irl-construction.png
	montage processing_speed_1-1.png processing_speed_1-2.png -tile 2x -gravity south -geometry +1+1 -bordercolor white -border 5x5 processing-speed.png

build:
	pandoc --standalone --filter ~/.cabal/bin/pandoc-crossref --template=rescience-template.tex --latex-engine=xelatex --biblatex --bibliography=article.bib -M "crossrefYaml=crossref.yaml" --output article.tex article.md
	xelatex article
	biber article
	xelatex article
	xelatex article
