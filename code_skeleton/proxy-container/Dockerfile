FROM python:3.9-slim

WORKDIR /app

RUN pip install pandas xgboost scikit-learn google-cloud-storage rdkit-pypi gcsfs "numpy<2"

COPY train_proxy.py .
COPY train_proxy_vertex.py .
COPY predict_proxy_vertex.py .
COPY score_library.py .

# Entrypoint can be overridden by the pipeline task
ENTRYPOINT ["python"]
