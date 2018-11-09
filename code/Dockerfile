FROM continuumio/anaconda3:latest

# Update the image since it sometimes are outdated
RUN conda update conda

# Install Neuron
ENV LANG=C.UTF-8

RUN apt-get update; apt-get install -y automake libtool build-essential openmpi-bin libopenmpi-dev \
                                       libncurses5-dev libreadline-dev libgsl0-dev cmake


ENV HOME=/home/docker
ENV VENV=$HOME/simulators
RUN mkdir $HOME; mkdir $HOME/packages; mkdir $VENV
ENV PATH=$PATH:$VENV/bin

ENV NRN_VER=7.6
ENV NRN=nrn-$NRN_VER

WORKDIR $HOME/packages
# Once neurons certificate is properly set up use add again
ADD http://www.neuron.yale.edu/ftp/neuron/versions/v$NRN_VER/$NRN.tar.gz .
RUN tar xzf $NRN.tar.gz; rm $NRN.tar.gz


RUN mkdir $VENV/build; mkdir $VENV/build/$NRN; mkdir $VENV/bin

WORKDIR $VENV/build/$NRN
RUN $HOME/packages/$NRN/configure --with-paranrn --with-nrnpython=python --disable-rx3d --without-iv --prefix=$VENV
RUN make
RUN make install
RUN cd src/nrnpython; python setup.py install
RUN cd $VENV/bin; ln -s ../x86_64/bin/nrnivmodl; ln -s ../x86_64/bin/nrngui; ln -s ../x86_64/bin/nrnoc; ln -s ../x86_64/bin/nrniv


WORKDIR /home/docker


# Install uncertainpy dependencies
RUN apt-get update --fix-missing
RUN apt-get -y install xvfb

RUN pip install xvfbwrapper

# Make sure matplotlib uses agg
RUN mkdir .config/
RUN mkdir .config/matplotlib
RUN echo "backend : Agg" >> .config/matplotlib/matplotlibrc


# Install latex, required for plotting with tex
RUN apt-get -y install texlive
RUN apt-get -y install texlive-latex-extra
RUN apt-get -y install dvipng


# install specific versions of packages used
RUN pip install --upgrade pip
RUN pip install numpy==1.15.2
RUN pip install matplotlib==3.0.0
RUN pip install uncertainpy==1.1.4
RUN pip install chaospy==2.3.5
RUN pip install tqdm==4.28.1


ENTRYPOINT cd code; nrnivmodl; bash


