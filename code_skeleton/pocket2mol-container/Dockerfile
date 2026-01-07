FROM continuumio/miniconda3

WORKDIR /app
COPY env_cuda113.yml .
RUN conda env create -f env_cuda113.yml
# Install system build tools (required for PyTorch extensions)
RUN apt-get update && apt-get install -y build-essential

RUN conda install -n Pocket2Mol -c conda-forge easydict pyyaml python-lmdb xgboost google-cloud-storage

# Explicitly install PyTorch first to satisfy dependencies (Match 1.10.0 for PyG wheels)
RUN /opt/conda/envs/Pocket2Mol/bin/pip install torch==1.10.0+cu113 -f https://download.pytorch.org/whl/cu113/torch_stable.html
# Then install PyG dependencies
RUN /opt/conda/envs/Pocket2Mol/bin/pip install torch-scatter torch-sparse torch-cluster torch-spline-conv -f https://data.pyg.org/whl/torch-1.10.0+cu113.html
# Finally install PyG
RUN /opt/conda/envs/Pocket2Mol/bin/pip install torch-geometric==2.0.4

# Copy the script and make it executable
COPY entrypoint.sh /usr/local/bin/
RUN chmod +x /usr/local/bin/entrypoint.sh

COPY batch_sample.sh .
COPY train.py .
COPY train_finetune.py .
COPY run_finetune.py .
COPY sample.py .
COPY sample_for_pdb.py .
COPY . .

# Set the entrypoint to our script.
# We DO NOT set CMD. We let GCP Batch provide the CMD.
# ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
