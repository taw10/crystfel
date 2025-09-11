FROM debian:trixie

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
   libgdk-pixbuf-2.0-dev \
   libfftw3-dev \
   git \
   flex \
   bison \
   libzmq3-dev \
   libmsgpack-dev \
   libcjson-dev \
   python3-dev \
   python3-pip \
   pandoc \
   unzip \
   wget \
   curl \
   ninja-build \
   gfortran \
   meson \
   && apt-get clean && rm -rf /var/lib/apt/lists/*

WORKDIR /home/crystfel-build

# Mosflm
RUN wget -nv https://www.mrc-lmb.cam.ac.uk/mosflm/mosflm/ver740/pre-built/mosflm-linux-64-noX11.zip
RUN unzip mosflm-linux-64-noX11.zip
RUN mv mosflm-linux-64-noX11 /usr/local/bin/mosflm

# CrystFEL
RUN git clone https://gitlab.desy.de/thomas.white/crystfel.git
RUN cd crystfel && meson setup build -Dprefix=/usr/local
RUN cd crystfel && meson compile -C build
RUN cd crystfel && meson test -C build
RUN cd crystfel && meson install -C build

## Stage 2
FROM debian:trixie-slim
RUN apt-get update && apt-get install -y \
   libhdf5-310 \
   libgsl28 \
   libgtk-3-0t64 \
   libcairo2 \
   libpango-1.0-0 \
   libgdk-pixbuf-2.0-0 \
   libfftw3-double3 \
   libpng16-16t64 \
   libzmq5 \
   libmsgpack-c2 \
   libcjson1 \
   && apt-get clean && rm -rf /var/lib/apt/lists/*
COPY --from=0 /usr/local /usr/local

# Environment variable needed for CrystFEL GUI, get_hkl and Mosflm
# (installed by Meson wrap for libccp4c)
ENV SYMINFO=/usr/local/share/ccp4/syminfo.lib

RUN ldconfig
