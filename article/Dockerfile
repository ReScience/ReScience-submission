FROM ivotron/texlive:20160320-1
MAINTAINER Ivo Jimenez <ivo.jimenez@gmail.com>

RUN apt-get -yq update && apt-get install -y python-pip wget
RUN wget https://github.com/jgm/pandoc/releases/download/1.19.2.1/pandoc-1.19.2.1-1-amd64.deb && \
    dpkg -i pandoc-1.19.2.1-1-amd64.deb && \
    wget https://github.com/lierdakil/pandoc-crossref/releases/download/v0.2.6.0/linux-ghc8-pandoc-1-19.tar.gz -q -O - | tar xz && \
    mv pandoc-crossref /usr/bin/ && \
    pip install pandocfilters && \
    apt-get clean -y && \
    rm -rf pandoc-1.19.2.1-1-amd64.deb /var/lib/apt/lists/* /tmp/* /var/tmp/*
RUN apt-get update
RUN apt-get install fonts-lmodern
RUN apt-get -y install texlive-fonts-extra
RUN apt-get -y install fonts-font-awesome
RUN apt-get update
RUN apt-get -y install biber --fix-missing
WORKDIR /root
ENTRYPOINT ["./build-in-docker.sh"]
#CMD ["--help"]
