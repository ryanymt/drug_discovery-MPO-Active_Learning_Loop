# Start from an official Python base image
FROM python:3.10-slim

# Set a working directory
WORKDIR /app

# Install system dependencies that some Python libraries might need
RUN apt-get update && apt-get install -y \
    build-essential \
    openbabel \
    && rm -rf /var/lib/apt/lists/*

# Install the necessary Python libraries
# MDAnalysis for trajectory analysis, Scikit-learn for clustering
RUN pip install --no-cache-dir \
    mdanalysis \
    scikit-learn \
    pandas \
    gcsfs

# Set a default entrypoint
#ENTRYPOINT ["/bin/bash"]
