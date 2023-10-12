FROM ubuntu:20.04
MAINTAINER Téo Lemane teo.lemane@inria.fr

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends && apt-get clean

RUN apt-get install -y git cmake gcc g++ libgsl-dev libopenblas-dev \
    liblapacke-dev libbz2-dev zlib1g zlib1g-dev python3

RUN cd /opt \
    && git clone --recursive https://github.com/tlemane/kmdiff \
    && cd kmdiff \
    && ./install.sh

RUN cd /opt/kmdiff && chmod +x ./kmdiff_build/bin/*

WORKDIR /tmp

ENTRYPOINT ["/opt/kmdiff/kmdiff_build/bin/kmdiff"]
