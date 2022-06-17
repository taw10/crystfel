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
   && apt-get clean && rm -rf /var/lib/apt/lists/*

RUN cd /usr/local/include && ln -sf eigen3/Eigen Eigen
RUN cd /usr/local/include && ln -sf eigen3/unsupported unsupported

WORKDIR /root
RUN pip3 install meson


RUN curl -L -O https://github.com/ninja-build/ninja/releases/download/v1.10.2/ninja-linux.zip
#wget https://github-releases.githubusercontent.com/ninja-linux.zip
RUN unzip ninja-linux.zip
RUN mv ninja /usr/local/bin

WORKDIR /home/crystfel-build
# Xgandalf
RUN git clone https://gitlab.desy.de/thomas.white/xgandalf.git
RUN mkdir xgandalf/build
RUN cd xgandalf && meson build/
RUN cd xgandalf && ninja -C build/
RUN cd xgandalf && ninja -C build/ install

# Fastdiffractionimageprocessing
RUN git clone https://stash.desy.de/scm/~gevorkov/fastdiffractionimageprocessing.git
RUN mkdir fastdiffractionimageprocessing/build
RUN cd fastdiffractionimageprocessing/build && cmake ..
RUN cd fastdiffractionimageprocessing && make -C build
RUN cd fastdiffractionimageprocessing && make -C build install

# Mosflm
RUN wget https://www.mrc-lmb.cam.ac.uk/mosflm/imosflm/ver740/downloads/imosflm-7.4.0-linux-64.zip
RUN unzip imosflm-7.4.0-linux-64.zip
RUN mv imosflm /usr/local/
RUN ln -sf ../imosflm/bin/mosflm /usr/local/bin/mosflm

# CrystFEL
RUN git clone --branch container https://gitlab.desy.de/silvan.schoen/crystfel.git
RUN mkdir crystfel/build
RUN cd crystfel && meson build/
RUN cd crystfel && ninja -C build/
RUN cd crystfel && ninja -C build/ install

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
ENV CLIBD=/usr/local/imosflm
ENV CINCL=/usr/local/imosflm
ENV CCP4_SCR=/usr/local/imosflm/src
RUN ldconfig

