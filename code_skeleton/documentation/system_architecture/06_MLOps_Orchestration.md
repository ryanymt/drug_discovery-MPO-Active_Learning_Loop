# MLOps & Orchestration

The platform employs a "Hybrid Orchestration" model, leveraging Google Cloud's managed services for heavy lifting while maintaining flexibility through Python-based control planes.

## 1. Compute Layer: Google Cloud Batch (The Workhorse)
*   **Role:** Massively Parallel Execution.
*   **Workloads:**
    *   **Generation:** 100k Pocket2Mol inference tasks.
    *   **Screening:** 70k RDKit/Gnina/TxGemma tasks.
    *   **Oracle:** 1000s of GROMACS MM-GBSA simulations.
*   **Configuration:** Defined in JSON specs (`workspace/*.json`) using the `gcloud batch` schema.
*   **Why Batch?** Traditional Kubernetes (GKE) is complex to tune for "burst" workloads. Batch provides ephemeral, "Serverless HPC" that scales to 0 when idle.

## 2. Training Layer: Vertex AI Custom Training
*   **Role:** Model Training & Fine-Tuning.
*   **Workloads:**
    *   **XGBoost:** Training the Proxy Model on Oracle data.
    *   **Pocket2Mol:** Fine-tuning the Generative Model on "Elite" candidates (A100 GPUs).
*   **Integration:** Jobs are submitted via `gcloud ai custom-jobs create`, ensuring full experiment tracking and metadata logging in Vertex AI.

## 3. Container Strategy (Artifact Registry)
We maintain a suite of specialized, versioned Docker images:

*   **`pocket2mol-retrain:v2`**: The Generative Engine.
    *   **Libraries**: PyTorch + PyG + RDKit.
    *   **Use Case**: Generation (Phase 1/10) & Fine-Tuning (Phase 9).
    
*   **`rdkit-gnina:latest`**: The Screening Workhorse.
    *   **Libraries**: RDKit, Gnina (Binary), OpenBabel, Scikit-Learn, XGBoost.
    *   **Use Case**: Docking, PhysChem Calculation, Proxy Inference (Phase 2/10b).

*   **`gromacs-mmpbsa:v1`**: The Oracle.
    *   **Libraries**: GROMACS 2023, AmberTools 22, Acpype.
    *   **Use Case**: MM-GBSA Binding Free Energy (Phase 4).

*   **`txgemma:v2`**: Large Language Model.
    *   **Libraries**: PyTorch, Transformers, vLLM.
    *   **Use Case**: Clinical Toxicity Classification (Phase 2/10b).

*   **`cloud-sdk:latest`**: Utility.
    *   **Libraries**: Google Cloud SDK (`gcloud`, `gsutil`, `bq`).
    *   **Use Case**: Data Prep, Consolidation, Manifest Generation.

## 4. The "Script Injection" Pattern
*   **Concept:** "Infrastructure as Code, Logic as Data".
*   **Problem:** Rebuilding 10GB Docker images for every minor script change is slow.
*   **Solution:**
    *   **Base Containers:** Hold heavy dependencies (CUDA, PyTorch, GROMACS).
    *   **Driver Scripts:** Python/Bash logic (`run_mmpbsa.sh`, `train_xgboost.py`) is stored in GCS (`gs://ryanymt/input/scripts/`).
    *   **Runtime:** Jobs download the latest script at startup.
    *   **Benefit:** We can patch logic in seconds without touching the infrastructure.

## 5. Metadata & Tracking
*   **Source of Truth:** The `molecule_registry` in BigQuery tracks every molecule's lifecyle.
*   **Run IDs:** Every batch job is tagged with a `run_id` (e.g., `prod_100k_v1`, `cycle2_10k`).
*   **Lineage:** We can trace a molecule from `global_id` (Generation) $\to$ `screening_results` (Filter) $\to$ `final_affinity` (Oracle) $\to$ `elite_candidates` (Training).
