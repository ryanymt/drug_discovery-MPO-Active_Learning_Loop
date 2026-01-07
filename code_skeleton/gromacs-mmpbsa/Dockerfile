# Start with the full NVIDIA CUDA development image for maximum compatibility
FROM nvidia/cuda:12.1.0-devel-ubuntu22.04

# Set non-interactive frontend to avoid prompts
ENV DEBIAN_FRONTEND=noninteractive

# 1. Install System Dependencies for GROMACS Build
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    wget \
    curl \
    tar \
    git \
    libfftw3-dev \
    dos2unix \
    && rm -rf /var/lib/apt/lists/*

# 2. Build GROMACS 2023.3 from Source (Ensures GPU/CUDA support)
WORKDIR /opt
RUN wget https://ftp.gromacs.org/gromacs/gromacs-2023.3.tar.gz && \
    tar xzvf gromacs-2023.3.tar.gz && \
    cd gromacs-2023.3 && \
    mkdir build && \
    cd build && \
    cmake .. -DGMX_BUILD_OWN_FFTW=ON \
    -DREGRESSIONTEST_DOWNLOAD=OFF \
    -DGMX_GPU=CUDA \
    -DCMAKE_INSTALL_PREFIX=/opt/gromacs \
    -DGMX_CUDA_TARGET_SM="75;80;86;89;90" && \
    make -j$(nproc) && \
    make install && \
    rm -rf /opt/gromacs-2023.3*

# 3. Install Micromamba for AmberTools (Solves HDF5/NetCDF ABI Hell)
# We install micromamba to /bin and set root prefix to /opt/conda
RUN curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xj -C /bin/ --strip-components=1 bin/micromamba

ENV MAMBA_ROOT_PREFIX=/opt/conda

# Create a robust conda environment
# We install acpype (which pulls ambertools) and other essentials
# Added mpi4py via conda to avoid compilation issues
RUN micromamba create -n fep_env \
    -c conda-forge -c bioconda \
    python=3.10 \
    acpype \
    ambertools \
    openbabel \
    numpy \
    mpi4py \
    pymol-open-source \
    && micromamba clean -a -y

# Install gmx_MMPBSA
# Note: mpi4py is now in conda.
RUN micromamba run -n fep_env pip install gmx_MMPBSA

# 4. Environment Setup
# Add GROMACS and Conda to PATH
ENV PATH="/opt/gromacs/bin:/opt/conda/envs/fep_env/bin:${PATH}"
# Set convenient aliases or variables if needed
ENV GMX_BIN="/opt/gromacs/bin/gmx"

# 5. Copy Execution Scripts
WORKDIR /opt
COPY run_gromacs.sh /usr/local/bin/run_gromacs.sh
COPY run_fep.sh /usr/local/bin/run_fep.sh
COPY fep_setup.py /usr/local/bin/fep_setup.py
COPY fep /opt/fep

# Fix permissions
RUN dos2unix /usr/local/bin/run_gromacs.sh /usr/local/bin/run_fep.sh && \
    chmod +x /usr/local/bin/run_gromacs.sh /usr/local/bin/run_fep.sh

# Default command
CMD ["/bin/bash"]