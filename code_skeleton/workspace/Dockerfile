# CORRECTED: Use the 'devel' image, which includes the full CUDA toolkit
# needed for compiling GPU-enabled applications from source.
FROM nvidia/cuda:12.1.0-devel-ubuntu22.04

# Avoid interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Install dependencies for GROMACS: build tools, git, and libraries
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    wget \
    git \
    libfftw3-dev \
    liblapack-dev \
    && rm -rf /var/lib/apt/lists/*

# Set a working directory
WORKDIR /opt

# Download, compile, and install GROMACS
# This process is unchanged, but will now succeed with the correct base image
RUN wget https://ftp.gromacs.org/gromacs/gromacs-2023.3.tar.gz && \
    tar xfz gromacs-2023.3.tar.gz && \
    cd gromacs-2023.3 && \
    mkdir build && \
    cd build && \
    cmake .. -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=ON -DGMX_GPU=CUDA -DCMAKE_INSTALL_PREFIX=/opt/gromacs && \
    make -j$(nproc) && \
    make install

# Add GROMACS to the system's PATH so it can be called directly
ENV PATH="/opt/gromacs/bin:${PATH}"

# Set a default entrypoint for the container
ENTRYPOINT ["/bin/bash"]
