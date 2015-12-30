
\[Re\] Least-cost modelling on irregular landscape graphs, J. Stachelek, Re**Science**, 2015.
  
**A reference implementation of** *Least-cost modelling on irregular landscape graphs*, T. Etherington, 27, 2012.

**Keywords**: Least-cost path, Delaunay triangulation, Graph Theory

##Prequisites 
 * git
 * R
 * GNU make
 * Tex
 * pandoc + pandoc-crossref

##PDF Build Instructions (tested on Ubuntu 14.04 and 15.10)
In a console, type:

```
$ sudo apt-get update
$ sudo apt-get -y install r-base
$ sudo apt-get install git
```

###devtools + dependencies
```
$ sudo apt-get install libcurl4-gnutls-dev
$ sudo apt-get install libxml2-dev
$ sudo apt-get install libssl-dev
$ sudo Rscript -e "install.packages('devtools', repos = 'https://cran.rstudio.com)"
```

###irlgraph + dependencies
```
$ sudo apt-get install libgeos-dev
$ sudo Rscript -e "devtools::install_github('jsta/irlgraph', dependencies = TRUE)"
```

###pandoc + dependencies
```
$ sudo Rscript -e "install.packages('rmarkdown', repos = 'https://cran.rstudio.com)"
$ sudo apt-get install texlive texlive-latex-extra texlive-xetex texlive-fonts-extra
$ sudo apt-get install texlive-bibtex-extra biber
$ sudo wget https://github.com/jgm/pandoc/releases/download/1.15.2/pandoc-1.15.2-1-amd64.deb
$ sudo dpkg -i pandoc-1.15.2-1-amd64.deb
$ sudo apt-get install cabal-install
$ cabal update
$ cabal install pandoc-crossref
```

###imagemagick
```
$ sudo apt-get install imagemagick
```

###checkout source and build
```
$ git clone https://github.com/jsta/ReScience-submission.git
$ cd ReScience-submission
$ git checkout STACHELEK
$ make all
```