# 01. Deployment Prerequisites

This document outlines the infrastructure and permissions required to deploy the "ChemOps" Active Learning Pipeline in a new GCP environment.

## 1. Google Cloud Project Setup
*   **APIs to Enable**:
    *   `aiplatform.googleapis.com` (Vertex AI)
    *   `batch.googleapis.com` (Cloud Batch)
    *   `compute.googleapis.com` (Compute Engine)
    *   `containerregistry.googleapis.com` (Container Registry) or `artifactregistry.googleapis.com`
    *   `storage-component.googleapis.com` (Cloud Storage)

## 2. Service Accounts & IAM
The pipeline requires a Service Account (SA) or user account to orchestrate jobs and accessing GCS.
*   **Recommended SA**: `vertex-ai-pipeline-sa@{project-id}.iam.gserviceaccount.com`
*   **Required Roles**:
    *   `Vertex AI User` (to submit pipelines)
    *   `Batch Job Editor` (to submit Cloud Batch jobs)
    *   `Storage Object Admin` (to read/write data)
    *   `Service Account User` (to attach itself to Compute instances)
    *   `Logs Writer` (for Cloud Logging)

## 3. Storage (GCS)
A central bucket is required for all inputs, outputs, and intermediate scripts.
*   **Bucket Name Constraint**: The pipeline code references `gs://drug-discovery-mvp-docking-results`.
    *   *Note*: If deploying to a new bucket, update `BASE_GCS_PATH` in `pipeline_e2e.py` and `job_spec` json templates.
*   **Secondary Bucket**: `gs://drug-discovery-mvp-input-data`
    *   **CRITICAL**: Used by `Pocket2Mol` to load the pretrained model (`models/pretrained.pt`) and target PDB (`4yhj.pdb`).
    *   *Action*: Ensure this bucket exists and contains the model artifacts.
*   **Directory Structure**:
    ```text
    gs://{bucket}/
    ├── data/
    │   └── scripts/           # Python scripts copied to containers at runtime
    │       ├── join_results.py
    │       ├── select_active_learning.py
    │       ├── train_xgboost.py
    │       ├── select_elite.py
    │       ├── redock_batch.py
    │       ├── finetune_launcher.py
    │       └── mock_fep_batch.py
    ├── pipeline_root/         # Vertex AI Metadata
    └── output/                # Dynamic output per loop_id
    ```

## 4. Compute & Quotas
Ensure the region (e.g., `us-central1`) has sufficient quota for the following machine families:
*   **Mock Pipeline**:
    *   `E2` CPUs (Spot): At least 20 vCPUs (for parallel screening).
*   **Prod Pipeline (FEP & Screening)**:
    *   `G2` (NVIDIA L4): Request `NVIDIA_L4_GPUS` quota.
        *   Used by: Gnina, TxGemma, GROMACS FEP.
        *   Count: At least 4-8 GPUs for parallel processing.
    *   `A2` (NVIDIA A100): Request `NVIDIA_A100_GPUS` quota.
        *   Used by: Pocket2Mol Fine-Tuning.
        *   Count: 1 GPU.

## 5. Container Images
Ensure the following images are pushed to your project's Artifact Registry:
*   `pocket2mol:latest` (Generation)
*   `gnina:latest` (Screening)
*   `txgemma:latest` (Screening)
*   `openbabel:latest` / `rdkit:latest` (Filtering/Redocking)
*   `gromacs:latest` (FEP Oracle)
*   `pocket2mol-retrain:v2` (Fine-Tuning)
*   `proxy-model:latest` (XGBoost/Selection/Launcher - General Python)
