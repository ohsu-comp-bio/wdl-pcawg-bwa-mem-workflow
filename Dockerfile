FROM ubuntu:14.04

LABEL version="2.6.7" \
      description="ICGC-TCGA-PanCaner BWA Workflow"

MAINTAINER Adam Struck <strucka@ohsu.edu>

USER root

ENV OPT /opt/icgc-tcga-pancan-deps
ENV PATH $OPT/bin:$PATH
ENV PERL5LIB $OPT/lib/perl5:$PERL5LIB

RUN apt-get update \
    && apt-get dist-upgrade -y --force-yes \
    && apt-get install software-properties-common -y

RUN apt-get -yqq update && \
    apt-get -yqq install libreadline6-dev\
                         build-essential \
                         g++ \
                         gfortran \
                         apt-utils \
                         autoconf \
                         pkg-config \
                         curl \
                         wget \
                         tar \
                         time \
                         zlib1g-dev \
                         liblz-dev \
                         libncurses5-dev \
                         libgd2-noxpm-dev \
                         cpanminus \
                         libwww-perl \
                         libxml-dom-perl \
                         libossp-uuid-perl \ 
                         libjson-perl \
                         libxml-libxml-perl \
                         libtry-tiny-perl \
                         libxml-xpath-perl \
                         libboost-dev \
                         libboost-iostreams-dev \
                         libpstreams-dev \
                         libglib2.0-dev \
                         samtools \ 
                         tabix

RUN apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN mkdir -p /tmp/downloads $OPT/bin $OPT/etc $OPT/lib $OPT/share
WORKDIR /tmp/downloads

RUN cpanm --mirror http://cpan.metacpan.org -l $OPT File::ShareDir File::ShareDir::Install Bio::Root::Version Const::Fast Graph && \
     rm -rf ~/.cpanm

# bwa - 0.7.12
RUN curl -ksSL -o tmp.tar.gz --retry 10 https://github.com/lh3/bwa/archive/0.7.12.tar.gz && \
    tar --strip-components 1 -zxf tmp.tar.gz && \
    make  && \
    cp bwa $OPT/bin/. && \
    rm -rf *

# biobambam2 - 2.0.35
RUN curl -ksSL -o tmp.tar.gz --retry 10 https://github.com/gt1/biobambam2/releases/download/2.0.35-release-20160330111451/biobambam2-2.0.35-release-20160330111451-x86_64-etch-linux-gnu.tar.gz && \
    tar --strip-components 1 -zxf tmp.tar.gz && \
    cp -r bin/* $OPT/bin/. && \
    cp -r etc/* $OPT/etc/. && \
    cp -r lib/* $OPT/lib/. && \
    cp -r share/* $OPT/share/. && \
    rm -rf *

# htslib - used multiple times later
RUN curl -ksSL -o tmp.tar.gz --retry 10 https://github.com/samtools/htslib/archive/1.2.1.tar.gz && \
    mkdir /tmp/downloads/htslib && \
    tar -C /tmp/downloads/htslib --strip-components 1 -zxf tmp.tar.gz && \
    make -C /tmp/downloads/htslib && \
    rm -f /tmp/downloads/tmp.tar.gz

ENV HTSLIB /tmp/downloads/htslib

# legacy samtools
RUN curl -ksSL -o tmp.tar.gz --retry 10 https://github.com/samtools/samtools/archive/0.1.20.tar.gz && \
    mkdir /tmp/downloads/samtools && \
    tar -C /tmp/downloads/samtools --strip-components 1 -zxf tmp.tar.gz && \
    perl -i -pe 's/^CFLAGS=\s*/CFLAGS=-fPIC / unless /\b-fPIC\b/' samtools/Makefile && \
    make -C samtools && \
    cp samtools/samtools $OPT/bin/. && \
    export SAMTOOLS=/tmp/downloads/samtools && \
    cpanm --mirror http://cpan.metacpan.org -l $OPT Bio::DB::Sam && \
    rm -rf /tmp/downloads/samtools /tmp/downloads/tmp.tar.gz ~/.cpanm

# bam_stats + PCAP build
RUN curl -ksSL -o tmp.tar.gz --retry 10 https://github.com/ICGC-TCGA-PanCancer/PCAP-core/archive/v1.13.1.tar.gz && \
    mkdir /tmp/downloads/PCAP && \
    tar -C /tmp/downloads/PCAP --strip-components 1 -zxf tmp.tar.gz && \
    make -C /tmp/downloads/PCAP/c && \
    cp /tmp/downloads/PCAP/bin/bam_stats $OPT/bin/. && \
    make -C /tmp/downloads/PCAP/c clean && \
    cd /tmp/downloads/PCAP && \
    cpanm --mirror http://cpan.metacpan.org -l $OPT . && \
    cd /tmp/downloads && \
    rm -rf /tmp/downloads/PCAP /tmp/downloads/tmp.tar.gz ~/.cpanm

# copy workflow specific scripts
COPY ./scripts/* $OPT/bin/

WORKDIR /home/

VOLUME /output/

CMD /bin/bash
