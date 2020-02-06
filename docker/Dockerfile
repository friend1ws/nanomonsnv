FROM ubuntu:18.04
MAINTAINER Kenichi Chiba <kchiba@hgc.jp> 

RUN apt-get update && apt-get install -y \
    git \
    wget \
    bzip2 \
    make \
    gcc \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    python3 \
    python3-pip

RUN wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && \
    tar jxvf samtools-1.9.tar.bz2 && \
    rm samtools-1.9.tar.bz2 && \
    cd samtools-1.9/htslib-1.9 && ./configure --enable-s3 && make && make install && \
    cd ../ && ./configure --without-curses && make && make install

RUN pip3 install --upgrade pip
RUN pip3 install --upgrade setuptools

ENV HTSLIB_LIBRARY_DIR /usr/local/lib
ENV HTSLIB_INCLUDE_DIR /usr/local/include

RUN pip3 install requests
RUN pip3 install cython==0.29.11
RUN pip3 install scipy==1.1.0
RUN pip3 install statistics==1.0.3.5
RUN pip3 install numpy==1.15.1
RUN pip3 install awscli boto3

RUN wget https://github.com/pysam-developers/pysam/archive/v0.15.2.tar.gz && \
    tar -zxvf v0.15.2.tar.gz && \
    rm v0.15.2.tar.gz && \
    cd pysam-0.15.2 && python3 setup.py build install

ENV LD_LIBRARY_PATH /usr/local/lib:$LD_LIBRARY_PATH

RUN wget https://github.com/ncc-ccat-gap/nanomonsnv/archive/v0.0.4.tar.gz && \
    tar -zxvf v0.0.4.tar.gz && \
    rm v0.0.4.tar.gz && \
    cd nanomonsnv-0.0.4 && \
    python3 setup.py build install

RUN apt-get install tabix
