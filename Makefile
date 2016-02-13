
#update_code:
#	rsync -av --exclude=.git --exclude='*.Rproj' --exclude='eth2012code' --exclude='.Rbuildignore' --exclude='.gitignore' --exclude='irlgraph_cache' --exclude='irlgraph_files' -- exclude='irlgraph.pdf' /home/jose/R/scripts/irlgraph code/

code/irlgraph/vignettes/irlgraph.pdf: code/irlgraph/vignettes/irlgraph.Rmd
	Rscript -e 'rmarkdown::render("code/irlgraph/vignettes/irlgraph.Rmd")'

move_images: code/irlgraph/vignettes/irlgraph.pdf
	rsync -av code/irlgraph/vignettes/irlgraph_files/figure-latex/ article

panel-resize_images: move_images
	cd article && make images

build_article: panel-resize_images
	cd article && make build

all: panel-resize_images
	@echo "figures built"

clean:
	rm -rf code/irlgraph/vignettes/irlgraph_cache
	rm -rf code/irlgraph/vignettes/irlgraph_files
	rm code/irlgraph/vignettes/irlgraph.pdf
	
