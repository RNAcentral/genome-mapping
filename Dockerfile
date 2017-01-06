# Create a reproducible installation of the initial version of mapping sequences
# to genomes.
#
# All local dependencies are installed manually to mirror the production setup
# where Docker or yum are not available.

FROM centos:6.6

RUN yum install -y \
    curl \
    gcc \
    git \
    libpng \
    openssl \
    openssl-devel \
    tar \
    unzip \
    zlib-devel

RUN mkdir /genome-mapping
RUN mkdir /genome-mapping/local
RUN mkdir /genome-mapping/bin
RUN mkdir /genome-mapping/app

ENV LOC /genome-mapping/local

# Install Python
RUN \
    cd $LOC && \
    curl -OL http://www.python.org/ftp/python/2.7.11/Python-2.7.11.tgz && \
    tar -zxvf Python-2.7.11.tgz && \
    cd Python-2.7.11 && \
    PREFIX=$LOC/python-2.7.11/ && \
    export LD_RUN_PATH=$PREFIX/lib && \
    ./configure --prefix=$PREFIX  --enable-shared && \
    make && \
    make install && \
    cd $LOC && \
    rm -Rf Python-2.7.11 && \
    rm Python-2.7.11.tgz

# Install virtualenv
RUN \
    cd $LOC && \
    curl -OL  https://pypi.python.org/packages/source/v/virtualenv/virtualenv-15.0.0.tar.gz && \
    tar -zxvf virtualenv-15.0.0.tar.gz && \
    cd virtualenv-15.0.0 && \
    $LOC/python-2.7.11/bin/python setup.py install && \
    cd $LOC && \
    rm -Rf virtualenv-15.0.0.tar.gz && \
    rm -Rf virtualenv-15.0.0

# Create RNAcentral virtual environment
RUN \
    cd $LOC && \
    mkdir virtualenvs && \
    cd virtualenvs && \
    $LOC/python-2.7.11/bin/virtualenv genome-mapping --python=$LOC/python-2.7.11/bin/python

# # Install BLAT
# RUN \
#   cd $LOC && \
#   curl -OL http://hgwdev.cse.ucsc.edu/~kent/src/blatSrc36.zip && \
#   unzip blatSrc36.zip && \
#   cd blatSrc && \
#   make MACHTYPE=x86_64 BINDIR=$LOC/bin

# Install python dependencies
ADD requirements.txt $LOC
RUN \
    source $LOC/virtualenvs/genome-mapping/bin/activate && \
    pip install -r $LOC/requirements.txt

WORKDIR /genome-mapping/app
COPY . /genome-mapping/app

CMD cd /genome-mapping/app && ls
CMD cd /genome-mapping/app && py.tests -v
