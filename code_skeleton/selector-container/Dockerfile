FROM python:3.9-slim

WORKDIR /app

RUN pip install pandas google-cloud-storage

COPY select_top.py .

ENTRYPOINT ["python", "select_top.py"]
