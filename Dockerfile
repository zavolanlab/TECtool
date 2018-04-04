# docker build -t tectool .
FROM ubuntu:16.04

RUN apt-get update && \
    apt-get install -y vim wget curl zlib1g-dev git software-properties-common libcurl4-gnutls-dev libxml2-dev libssl-dev apt-transport-https libmariadb-client-lgpl-dev && \
    apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 && \
    add-apt-repository 'deb [arch=amd64,i386] https://cran.rstudio.com/bin/linux/ubuntu xenial/' && \
    apt-get update 

# install R
RUN apt-get install -y r-base && \
    Rscript -e 'install.packages(c("devtools", "optparse"), repos = "http://cran.us.r-project.org")' && \
    Rscript -e 'source("https://bioconductor.org/biocLite.R"); biocLite(c("rtracklayer", "Gviz", "biomaRt", "GenomicFeatures"), ask=FALSE);'

# install python 
RUN apt-get install -y build-essential python3 python3-dev python3-pip && \
    python3 -m pip install pip --upgrade && \
    python3 -m pip install wheel && \
    unlink /usr/bin/python && \
    ln -s /usr/bin/python3 /usr/bin/python

# fix vim
RUN cd $HOME && \
    git clone https://github.com/fgypas/.vim.git && \
    cp $HOME/.vim/vimrc /etc/vim/vimrc.local

# install bedtools
RUN cd $HOME && \
    wget https://github.com/arq5x/bedtools2/archive/v2.26.0.tar.gz && \
    tar xzvf v2.26.0.tar.gz && \
    rm v2.26.0.tar.gz && \
    cd bedtools2-2.26.0 && \
    make && \
    make install

# install tectool
RUN cd $HOME && \
    git clone https://git.scicore.unibas.ch/zavolan_public/TECtool.git && \
    cd TECtool && \
    pip install -r requirements.txt && \
    python3 setup.py install

## download test data
#RUN cd $HOME/TECtool && \
#    wget http://tectool.unibas.ch/data/test_data.tar.gz && \
#    tar xzvf test_data.tar.gz && \
#    rm test_data.tar.gz && \
#    cd test_data
