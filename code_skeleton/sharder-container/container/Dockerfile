
# Use a lightweight Python base image
FROM python:3.9-slim

# Set the working directory
WORKDIR /app

# Copy the sharder script into the container
COPY sharder.py .

# The script will be executed via command-line arguments in the Batch job
ENTRYPOINT ["python", "sharder.py"]
