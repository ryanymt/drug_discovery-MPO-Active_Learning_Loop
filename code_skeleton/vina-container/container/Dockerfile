# Start from a minimal Debian base image
FROM debian:bullseye-slim

# Install wget and dependencies
RUN apt-get update && apt-get install -y \
    wget \
    bzip2 \
    && rm -rf /var/lib/apt/lists/*

# Set a working directory
WORKDIR /opt

# Download and extract AutoDock Vina
RUN wget https://vina.scripps.edu/wp-content/uploads/sites/55/2020/12/autodock_vina_1_1_2_linux_x86.tgz && \
    tar xzvf autodock_vina_1_1_2_linux_x86.tgz

# Add the Vina binaries to the system's PATH
ENV PATH="/opt/autodock_vina_1_1_2_linux_x86/bin:${PATH}"

# Set a default entrypoint
#ENTRYPOINT ["/bin/bash"]