FROM python:3.9-slim

WORKDIR /app

RUN pip install google-cloud-aiplatform google-cloud-storage pandas

COPY suggest.py .

ENTRYPOINT ["python", "suggest.py"]
