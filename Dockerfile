FROM ubuntu:20.04

SHELL ["/bin/bash", "-c"]

ARG ENABLE_PBAM=ON
ARG ENABLE_PBSAM=ON
ARG ENABLE_PBAM_SPHINX=OFF
ARG ENABLE_PB_TESTING=OFF

RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y \
        build-essential \
        cmake \
        libopenblas-dev \
        python3-dev \
        cython \
        && \
    mkdir /src && \
    /bin/true

COPY cmake /src/cmake
COPY pb_shared /src/pb_shared
COPY pb_wrap /src/pb_wrap
COPY pbam /src/pbam
COPY pbsam /src/pbsam
COPY scripts /src/scripts
COPY CMakeLists.txt /src/CMakeLists.txt

RUN cd /src && \
    mkdir build && \
    cd build && \
    cmake -DCMAKE_INSTALL_PREFIX=/src/build -DENABLE_PBAM=${ENABLE_PBAM} -DENABLE_PBSAM=${ENABLE_PBSAM} -DENABLE_PBAM_SPHINX=${ENABLE_PBAM_SPHINX} -DENABLE_PB_TESTING=${ENABLE_PB_TESTING} .. && \
    make install && \
    /bin/true
