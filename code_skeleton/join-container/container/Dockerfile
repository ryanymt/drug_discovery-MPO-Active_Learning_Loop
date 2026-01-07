FROM python:3.9-slim

# Install system dependencies
RUN apt-get update && apt-get install -y \
    git \
    && rm -rf /var/lib/apt/lists/*

# Install Python libraries
RUN pip install --no-cache-dir \
    pandas \
    google-cloud-storage \
    google-cloud-bigquery \
    pyarrow \
    db-dtypes

WORKDIR /app

COPY join_results.py .

ENTRYPOINT ["python", "join_results.py"]

