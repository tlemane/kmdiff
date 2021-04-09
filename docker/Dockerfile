FROM ubuntu:latest
MAINTAINER Téo Lemane teo.lemane@inria.fr

ENV DEBIAN_FRONTEND=noninteractive 

RUN apt-get update && apt-get -y dist-upgrade \
    && apt-get install -y --no-install-recommends && apt-get clean

RUN apt-get install -y git cmake gcc g++ zlib1g zlib1g-dev python3

RUN cd /opt \
    && git clone --recursive https://github.com/tlemane/kmdiff \
    && cd kmdiff \
    && mkdir build \
    && cd build \
    && cmake .. \
    && make

RUN cd /opt/kmdiff && chmod +x ./bin/*

WORKDIR /tmp

ENTRYPOINT ["/opt/kmdiff/bin/kmdiff"]