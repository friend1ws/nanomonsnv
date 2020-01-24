FROM python:3.7.5-buster
MAINTAINER kchiba <kchiba@hgc.jp>

RUN apt-get -y update && \
    apt-get install -y wget unzip && \
    \
    wget -nc https://github.com/ncc-ccat-gap/nanomonsnv/archive/v0.0.1.zip && \
    unzip v0.0.1.zip && \
    rm -f v0.0.1.zip && \
    cd nanomonsnv-0.0.1 && \
    python setup.py build install && \
    \
    mkdir -p /tools && \
    \
    cd /tools && \
    echo "build samtools-1.9" && \
    wget -nc https://sourceforge.net/projects/samtools/files/samtools/1.9/samtools-1.9.tar.bz2 && \
    tar xjvf samtools-1.9.tar.bz2 && \
    rm -f samtools-1.9.tar.bz2 && \
    cd samtools-1.9 && \
    make && \
    make install && \
    \
    cd /tools && \
    echo "build htslib-1.9" && \
    wget -nc https://sourceforge.net/projects/samtools/files/samtools/1.9/htslib-1.9.tar.bz2 && \
    tar xvf htslib-1.9.tar.bz2 && \
    rm -f htslib-1.9.tar.bz2 && \
    cd htslib-1.9 && \
    make && \
    make install

RUN pip install --upgrade pip
RUN pip install scipy==1.1.0
RUN pip install statistics==1.0.3.5
RUN pip install pysam==0.15.2
RUN pip install numpy==1.15.1

CMD ["/bin/bash"]
