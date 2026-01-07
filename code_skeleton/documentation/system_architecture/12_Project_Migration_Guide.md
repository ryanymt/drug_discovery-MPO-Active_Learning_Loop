# Project Migration & "Warm-Up" Guide

**Objective:** Migrate all necessary Checkpoints, Data, Scripts, and Container Images from the current POC project (`lifescience-project-469915`) to a new Production Project.

**Prerequisites:**
*   You have access to both projects.
*   You have `gcloud`, `gsutil`, and `docker` installed locally.
*   You have roughly **50GB** of local disk space (for Docker images and checkpoint files).

---

## Part 1: GCS Asset Migration (Data & Scripts)

This section downloads the "Source of Truth" data to your local machine and prepares it for upload.

### 1.1 Download Assets Locally
Run this in a clean shell.

```bash
# 1. Create a local staging directory
mkdir -p ./migration_stage/input
mkdir -p ./migration_stage/scripts
mkdir -p ./migration_stage/models

# 2. Authenticate to SOURCE Project (Current)
gcloud config set project lifescience-project-469915
gcloud auth login

# 3. Download Protein Structure (Receptor)
gsutil cp gs://drug-discovery-mvp-input-data/4yhj.pdb ./migration_stage/input/
gsutil cp gs://drug-discovery-mvp-input-data/1aq1.pdb ./migration_stage/input/
gsutil cp gs://drug-discovery-mvp-input-data/1aki.pdb ./migration_stage/input/

# 4. Download Pre-trained Models
gsutil cp gs://drug-discovery-mvp-input-data/models/pretrained.pt ./migration_stage/models/

# 5. Download Pipeline Scripts
# (Note: Some are in 'input-data' and some in 'docking-results')
gsutil cp gs://drug-discovery-mvp-docking-results/data/scripts/*.py ./migration_stage/scripts/
gsutil cp gs://drug-discovery-mvp-input-data/scripts/predict_batched.py ./migration_stage/scripts/
```

### 1.2 Upload Assets to DESTINATION Project
Switch to your new project and upload.

```bash
# 1. Authenticate to NEW Project
export NEW_PROJECT_ID="YOUR_NEW_PROJECT_ID"
gcloud config set project $NEW_PROJECT_ID

# 2. Create Buckets
gsutil mb -p $NEW_PROJECT_ID -l us-central1 gs://${NEW_PROJECT_ID}-input
gsutil mb -p $NEW_PROJECT_ID -l us-central1 gs://${NEW_PROJECT_ID}-results

# 3. Upload Data
gsutil cp ./migration_stage/input/*.pdb gs://${NEW_PROJECT_ID}-input/
gsutil cp ./migration_stage/models/pretrained.pt gs://${NEW_PROJECT_ID}-input/models/

# 4. Upload Scripts
# Note: Pipeline expects scripts in a specific path "data/scripts"
gsutil cp ./migration_stage/scripts/*.py gs://${NEW_PROJECT_ID}-results/data/scripts/
```

---

## Part 2: Container Image Migration

We will pull the **exact verified versions** from the old Artifact Registry, re-tag them, and push to the new one.

### 2.1 Setup & Auth
```bash
# Authenticate Docker to gcr.io and pkg.dev
gcloud auth configure-docker us-central1-docker.pkg.dev
```

### 2.2 Pull Images (From Source)
```bash
export SOURCE_REG="us-central1-docker.pkg.dev/lifescience-project-469915/drug-discovery-containers"

# Pull the specific tags used in Production Runbook
docker pull ${SOURCE_REG}/pocket2mol:v3
docker pull ${SOURCE_REG}/pocket2mol-retrain:v2
docker pull ${SOURCE_REG}/gnina:latest
docker pull ${SOURCE_REG}/txgemma:v2
docker pull ${SOURCE_REG}/gromacs-mmpbsa:latest
docker pull ${SOURCE_REG}/proxy-model:latest
docker pull ${SOURCE_REG}/rdkit:latest
```

### 2.3 Retag & Push (To Destination)
```bash
export NEW_REG="us-central1-docker.pkg.dev/${NEW_PROJECT_ID}/drug-discovery-containers"

# 1. Create Repository in New Project
gcloud artifacts repositories create drug-discovery-containers \
    --repository-format=docker \
    --location=us-central1 \
    --description="Drug Discovery Pipeline Containers"

# 2. Retag and Push Loop
# Pocket2Mol (Generation)
docker tag ${SOURCE_REG}/pocket2mol:v3 ${NEW_REG}/pocket2mol:v3
docker push ${NEW_REG}/pocket2mol:v3

# Pocket2Mol (Retraining)
docker tag ${SOURCE_REG}/pocket2mol-retrain:v2 ${NEW_REG}/pocket2mol-retrain:v2
docker push ${NEW_REG}/pocket2mol-retrain:v2

# Gnina
docker tag ${SOURCE_REG}/gnina:latest ${NEW_REG}/gnina:latest
docker push ${NEW_REG}/gnina:latest

# TxGemma
docker tag ${SOURCE_REG}/txgemma:v2 ${NEW_REG}/txgemma:v2
docker push ${NEW_REG}/txgemma:v2

# Gromacs MM-GBSA
docker tag ${SOURCE_REG}/gromacs-mmpbsa:latest ${NEW_REG}/gromacs-mmpbsa:latest
docker push ${NEW_REG}/gromacs-mmpbsa:latest

# Proxy Model (XGBoost)
docker tag ${SOURCE_REG}/proxy-model:latest ${NEW_REG}/proxy-model:latest
docker push ${NEW_REG}/proxy-model:latest

# RDKit (Helper)
docker tag ${SOURCE_REG}/rdkit:latest ${NEW_REG}/rdkit:latest
docker push ${NEW_REG}/rdkit:latest
```

---

## Part 3: Pipeline Config Update

After uploading, you must update the JSON config files in your local `workspace/` to point to the new project and bucket names before running `pipeline_e2e.py`.

1.  **Replace Project ID**: Find/Replace `lifescience-project-469915` -> `$NEW_PROJECT_ID`.
2.  **Replace Buckets**:
    *   `drug-discovery-mvp-input-data` -> `${NEW_PROJECT_ID}-input`
    *   `drug-discovery-mvp-docking-results` -> `${NEW_PROJECT_ID}-results`

Run the pipeline compilation again:
```bash
python3 workspace/pipeline_e2e.py
```
