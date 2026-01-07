# Production Test Plan: 100k Molecule Active Learning Loop

**Objective:** Validate the full-scale Active Learning loop (100,000 molecules) using manual job submission before automating with Kubeflow.

**Prerequisites:**
*   GCP Project with **160+ L4 GPU Quota**.
*   GCS Buckets created: `gs://[PROJECT_ID]-results`, `gs://[PROJECT_ID]-input`.
*   BigQuery Dataset created: `bioops_platform`.
*   Docker Images pushed to Artifact Registry.

---

## Phase 1: Massive Generation
**Goal:** Generate 100,000 raw candidates.

### Job 1: Pocket2Mol Generation (Cloud Batch)
*   **Config File:** `jobs/01_generation_batch.json`
*   **Container:** `us-central1-docker.pkg.dev/[PROJECT_ID]/drug-discovery-containers/pocket2mol:v3`
*   **Task Count:** 100 (Each task generates 1,000 mols).
*   **Parallelism:** 50 (or Max Quota).
*   **Input:**
    *   `gs://[PROJECT_ID]-input/4yhj.pdb` (Protein structure).
    *   `gs://[PROJECT_ID]-input/models/pretrained.pt` (Checkpoint).
*   **Command:**
    ```bash
    # Note: Requires runtime config setup and env vars
    export LD_LIBRARY_PATH=/usr/local/nvidia/lib64:/opt/conda/envs/Pocket2Mol/lib:$LD_LIBRARY_PATH && \
    cp configs/sample_for_pdb.yml /tmp/run_config.yml && \
    sed -i 's|./ckpt/pretrained_Pocket2Mol.pt|/mnt/disks/models/models/pretrained.pt|g' /tmp/run_config.yml && \
    python sample_for_pdb.py \
      --config /tmp/run_config.yml \
      --pdb_path /mnt/disks/inputs/4yhj.pdb \
      --center ' 0.517,27.062,8.972' \
      --outdir /mnt/disks/outputs/cycle_01/raw_sdfs/shard_${BATCH_TASK_INDEX}
    ```
*   **Expected Output:** 100,000 `.sdf` files in `gs://[PROJECT_ID]-results/cycle_01/raw_sdfs/shard_*/`.

---

## Phase 2: Parallel Scoring (The Filter)
**Goal:** Score all 100,000 molecules for Affinity and QED.

### Job 2: RDKit + Gnina Scoring (Cloud Batch)
*   **Config File:** `workspace/gnina-baseline-batch.json`
*   **Container:** `us-central1-docker.pkg.dev/[PROJECT_ID]/drug-discovery-containers/gnina:latest`
*   **Task Count:** 20 (Each task processes 5,000 mols).
*   **Input:** `gs://[PROJECT_ID]-results/cycle_01/raw_sdfs/` (The output from Step 1).
*   **Command:**
    ```bash
    # Uses wrapper script that iterates through shards
    export LD_LIBRARY_PATH=/usr/local/nvidia/lib64:/usr/local/cuda/lib64:$LD_LIBRARY_PATH && \
    /bin/bash /mnt/disks/gcs/scripts/run_gnina_baseline.sh
    ```
*   **Expected Output:** 20 CSV files in `gs://[PROJECT_ID]-results/cycle_01/scores/`.

### Job 2b: TxGemma ADMET Scoring (Cloud Batch)
### Job 2b: TxGemma ADMET Scoring (Cloud Batch)
*   **Config File:** `workspace/txgemma-baseline-batch.json`
*   **Container:** `us-central1-docker.pkg.dev/[PROJECT_ID]/drug-discovery-containers/txgemma:v2`
*   **Task Count:** 20 (Same sharding).
*   **Parallelism:** 10 (GPU dependent).
*   **Input:** `gs://[PROJECT_ID]-results/cycle_01/raw_sdfs/` (or SMILES lists).
*   **Command:**
    ```bash
    # Runtime dependency installation required
    echo 'Starting Job...' && \
    echo 'Installing Transformers...' && \
    pip install transformers==4.42.4 --quiet --no-cache-dir && \
    echo 'Setup Complete. Copying script...' && \
    cp /mnt/disks/gcs/data/scripts/predict_batched.py /tmp/predict.py && \
    echo 'Running Inference...' && \
    python3 -u /tmp/predict.py \
      --input-file /mnt/disks/gcs/input/baseline/baseline_smiles.txt \
      --output-file /mnt/disks/gcs/input/baseline/analysis/txgemma_results_baseline.csv
    ```
*   **Expected Output:** 20 CSV files with `toxicity_score`.

---

## Phase 3: Data Ingestion & Selection
**Goal:** Centralize data and pick the "Oracle Set".

### Job 3: Joiner (Local Script or Small Job)
*   **Action:**
    1.  Download all 20 CSVs.
    2.  Concatenate into `cycle_01_all_scores.csv`.
    3.  Generate `molecule_hash` (SHA256 of SMILES).
    4.  Upload to BigQuery `screening_results`.
*   **Selection Logic (SQL):**
    ```sql
    SELECT * FROM `bioops_platform.screening_results`
    ORDER BY (0.7 * gnina_affinity + 0.3 * qed_score) DESC
    LIMIT 900
    UNION ALL
    SELECT * FROM `bioops_platform.screening_results` ORDER BY RAND() LIMIT 100
    ```
*   **Expected Output:** `oracle_candidates.csv` (1,000 rows) uploaded to GCS.

---

## Phase 4: The Oracle (MM-GBSA)
**Goal:** Get physics-based Ground Truth for the top 1,000.

### Job 4: GROMACS MM-GBSA (Cloud Batch)
### Job 4: GROMACS MM-GBSA (Cloud Batch)
*   **Config File:** `workspace/mmpbsa-test-job.json`
*   **Container:** `us-central1-docker.pkg.dev/[PROJECT_ID]/drug-discovery-containers/gromacs-mmpbsa:latest`
*   **Task Count:** 10 (Each task runs 100 simulations).
*   **Parallelism:** 10 (10 GPUs).
*   **Input:**
    *   `gs://[PROJECT_ID]-input/oracle_candidates.csv` (The 1k list).
    *   `gs://[PROJECT_ID]-input/4yhj.pdb` (Receptor).
*   **Command:**
    ```bash
    # Hotfix scripts at runtime and execute
    cp /mnt/disks/gcs/input/fep_setup.py /usr/local/bin/fep_setup.py && \
    cp /mnt/disks/gcs/input/run_mmpbsa.sh /usr/local/bin/run_mmpbsa.sh && \
    chmod +x /usr/local/bin/run_mmpbsa.sh /usr/local/bin/fep_setup.py && \
    export BATCH_TASK_INDEX=${BATCH_TASK_INDEX} && \
    /usr/local/bin/run_mmpbsa.sh --run --id batch_${BATCH_TASK_INDEX}
    ```
*   **Expected Output:** `gs://[PROJECT_ID]-results/cycle_01/oracle/mmpbsa_results.csv` (Aggregated or 10 parts).

---

## Phase 5: The "Brain" (XGBoost Training)
**Goal:** Learn the mapping from Structure -> MM-GBSA Energy.

### Job 5: Train Proxy (Vertex AI Custom Job)
*   **Config File:** `jobs/05_train_config.yaml`
*   **Container:** `us-central1-docker.pkg.dev/[PROJECT_ID]/drug-discovery-containers/proxy-model:latest`
*   **Machine:** `n1-standard-4` (No GPU needed for XGBoost on small data).
*   **Input:** `gs://[PROJECT_ID]-results/cycle_01/oracle/mmpbsa_results.csv` (The 1k labeled).
*   **Command:**
    ```bash
    python3 train_xgboost.py \
      --train /gcs/[BUCKET]/cycle_01/oracle/mmpbsa_results.csv \
      --output_model /gcs/[BUCKET]/models/cycle_01_xgboost.json
    ```
*   **Expected Output:** `gs://[PROJECT_ID]-results/models/cycle_01_xgboost.json`.

---

## Phase 6: Inference & Elite Selection
**Goal:** Score the 99,000 skipped molecules.

### Job 6: Batch Inference (Vertex AI or Batch)
*   **Container:** `proxy-model:latest`.
*   **Input:**
    *   `cycle_01_all_scores.csv` (The full 100k list).
    *   `cycle_01_xgboost.json`.
*   **Output:** `predicted_scores.csv`.
*   **Selection:** Filter for Top 10,000 based on *Predicted* MM-GBSA score. Save as `elite_5k.csv`.

---

## Phase 7: Fine-Tuning (The Loop Closure)
**Goal:** Update Pocket2Mol weights.

### Job 7: Retraining (Vertex AI Custom Job)
*   **Config File:** `jobs/07_finetune_config.yaml`
*   **Container:** `pocket2mol-retrain:v2`
*   **Machine:** `a2-highgpu-1g` (A100) or `n1-standard-4` (T4).
*   **Input:**
    *   `elite_5k.csv` (SMILES list).
    *   **Crucial Step:** The job script must first **Redock** these 5k SMILES to generate 3D labels (SDFs) *before* starting training.
*   **Command:**
    ```bash
    python3 run_finetune_v7.py \
      --elite_csv /gcs/[BUCKET]/elite_5k.csv \
      --base_ckpt /gcs/[BUCKET]/models/pretrained.pt \
      --output_ckpt /gcs/[BUCKET]/models/cycle_02.pt
    ```
*   **Expected Output:** `gs://[PROJECT_ID]-results/models/cycle_02.pt`.
