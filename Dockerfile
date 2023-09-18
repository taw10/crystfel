FROM debian:buster

RUN apt-get update && apt-get install -y \
   pkg-config \
   cmake \
   build-essential \
   libhdf5-dev \
   libgsl-dev \
   libgtk-3-dev \
   libcairo2-dev \
   libeigen3-dev \
   libpango1.0-dev \
   libgdk-pixbuf2.0-dev \
   libfftw3-dev \
   libncurses-dev \
   libpng-dev \
   libtiff5-dev \
   git \
   flex \
   bison \
   libzmq3-dev \
   libmsgpack-dev \
   python3-dev \
   python3-pip \
   unzip \
   wget \
   curl \
   ninja-build \
   gfortran \
   && apt-get clean && rm -rf /var/lib/apt/lists/*

WORKDIR /root
RUN pip3 install meson

WORKDIR /home/crystfel-build

# Mosflm
RUN wget -nv https://www.mrc-lmb.cam.ac.uk/mosflm/mosflm/ver740/pre-built/mosflm-linux-64-noX11.zip
RUN unzip mosflm-linux-64-noX11.zip
RUN mv mosflm-linux-64-noX11 /usr/local/bin/mosflm

# CrystFEL
RUN git clone https://gitlab.desy.de/thomas.white/crystfel.git
RUN cd crystfel && meson build -Dprefix=/usr/local
RUN cd crystfel && ninja -C build
RUN cd crystfel && ninja -C build test
RUN cd crystfel && ninja -C build install

## Stage 2
FROM debian:buster-slim
RUN apt-get update && apt-get install -y \
   libhdf5-103 \
   libgsl23 \
   libgtk-3-0 \
   libcairo2 \
   libpango1.0 \
   libgdk-pixbuf2.0 \
   libfftw3-double3 \
   libncurses6 \
   libpng16-16 \
   libtiff5 \
   libzmq5 \
   libmsgpackc2 \
   && apt-get clean && rm -rf /var/lib/apt/lists/*
COPY --from=0 /usr/local /usr/local

# Environment variable needed for CrystFEL GUI and Mosflm
# The file is installed by libccp4c, a wrapped subproject of CrystFEL
ENV SYMINFO=/usr/share/ccp4/syminfo.lib

RUN ldconfig
