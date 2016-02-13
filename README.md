
\[Re\] Least-cost modelling on irregular landscape graphs, J. Stachelek, Re**Science**, 2015.
  
**A reference implementation of** *Least-cost modelling on irregular landscape graphs*, T. Etherington, 27, 2012.

**Keywords**: Least-cost path, Delaunay triangulation, Graph Theory

##Figure Build Instructions (tested on Ubuntu 14.04 and 15.10)
###Prequisites 
 * git
 * R
 * GNU make
 * pandoc

In a console, type:

```
$ echo deb https://cran.rstudio.com/bin/linux/ubuntu trusty/ >> /etc/apt/sources.list
$ sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9  
$ sudo apt-get update
$ sudo apt-get -y install r-base
$ sudo apt-get install git
```

###devtools + dependencies

```
$ sudo apt-get install libcurl4-gnutls-dev
$ sudo apt-get install libxml2-dev
$ sudo apt-get install libssl-dev
$ sudo Rscript -e "install.packages('devtools', repos = 'https://cran.rstudio.com')"
```

###irlgraph + dependencies
```
$ sudo apt-get install libgeos-dev
$ sudo Rscript -e "devtools::install_github('jsta/irlgraph', dependencies = TRUE)"
```

###pandoc + dependencies
```
$ sudo Rscript -e "install.packages('rmarkdown', repos = 'https://cran.rstudio.com')"
$ sudo wget https://github.com/jgm/pandoc/releases/download/1.15.2/pandoc-1.15.2-1-amd64.deb
$ sudo dpkg -i pandoc-1.15.2-1-amd64.deb
$ sudo apt-get install texlive texlive-latex-extra
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

##Docker Build Instructions
###Prequisites 
 * docker

In a console, type:

```
$ docker pull jsta/irlgraph-test
$ docker cp $(docker create jsta/irlgraph-test):ReScience-submission/article/article.pdf .
```

##PDF Build Instructions

```
$ sudo apt-get install texlive-xetex texlive-fonts-extra
$ sudo apt-get install texlive-bibtex-extra biber
$ sudo apt-get install cabal-install
$ cabal update
$ cabal install pandoc-crossref
$ make all
$ make build_article
```