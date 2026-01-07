# Use a base image with CUDA 12.1 and cuDNN 8 for GPU support
FROM nvidia/cuda:12.1.0-cudnn8-devel-ubuntu22.04

# Avoid prompts from apt
ENV DEBIAN_FRONTEND=noninteractive

# Install Python and pip
RUN apt-get update && apt-get install -y \
    python3.9 \
    python3-pip \
    && rm -rf /var/lib/apt/lists/*

# Install required Python libraries
RUN pip3 install --no-cache-dir \
    torch==2.1.0 \
    transformers==4.35.2 \
    accelerate==0.25.0 \
    pandas==2.1.3 \
    'google-cloud-storage' \
    safetensors \
    bitsandbytes==0.41.3 \
    scipy \
    huggingface_hub

# Set the working directory
WORKDIR /app

# Copy the prediction script into the container
COPY predict.py .

# The script will be executed via command-line arguments in the Batch job
ENTRYPOINT ["python3", "predict.py"]
