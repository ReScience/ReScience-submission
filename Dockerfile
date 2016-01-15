#not installed texlive-fonts-extra 

FROM achubaty/r-spatial-base

MAINTAINER Joseph Stachelek <jstachel@sfwmd.gov>

RUN apt-get update \
	&& apt-get install -y \
		git \
		rsync \
		imagemagick \
		texlive \
		texlive-latex-extra \
		texlive-xetex \
		texlive-bibtex-extra \
		fonts-roboto-hinted \
		fonts-roboto \
		texlive-fonts-extra \
		biber \
		cabal-install

RUN Rscript -e 'devtools::install_github("jsta/irlgraph", dependencies = TRUE)'
RUN Rscript -e "install.packages('rmarkdown', repos = 'https://cran.rstudio.com')"

RUN wget https://github.com/jgm/pandoc/releases/download/1.15.2/pandoc-1.15.2-1-amd64.deb \
	&& dpkg -i pandoc-1.15.2-1-amd64.deb

RUN cabal update \
	&& cabal install pandoc-crossref-0.1.5.6

RUN git clone https://github.com/jsta/ReScience-submission.git \
	&& cd ReScience-submission \
		&& git checkout STACHELEK \
			&& make all



