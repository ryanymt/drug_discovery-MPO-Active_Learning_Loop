FROM python:3.9-slim

RUN pip install pandas google-cloud-bigquery google-cloud-storage

WORKDIR /app
COPY join.py /app/join.py

ENTRYPOINT ["python", "/app/join.py"]
