# 02. Pipeline Overview & Architecture

## Pipeline Versions
We maintain two versions of the end-to-end pipeline, generated from the single source of truth `pipeline_e2e.py`.

### 1. Mock Pipeline (`00-mock-pipeline.json`)
*   **Purpose**: Rapid validation of the orchestration logic, data flow, and error handling.
*   **Differences**: Unlocks the "Oracle" step to run `mock_fep_batch.py` (simulated scores) instead of GROMACS. Uses lightweight CPU instances for most steps.
*   **Validation Command**:
    ```bash
    python3 workspace/submit_pipeline.py # Submits '00-mock-pipeline.json'
    ```

### 2. Prod Pipeline (`01-prod-pipeline.json`)
*   **Purpose**: Real scientific production run.
*   **Differences**: Uses `01-prod-fep-task-spec.json` (Real GROMACS) and GPU-accelerated screening.
*   **Note**: Ensure quota availability before running.

## Active Learning Flow (Loop Architecture)
The pipeline executes a "Closed Loop" cycle consisting of 9 steps:

1.  **Generation** (`pocket2mol`): Generates SMILES strings.
2.  **Filtering** (Parallel):
    *   `Gnina` (Docking)
    *   `TxGemma` (Target Affinity/Toxicity)
    *   `RDKit` (QED, LogP, SA)
3.  **Data Joiner**: Merges all filter results into `joined_results.csv`.
4.  **Selection (Active Learning)**:
    *   Splits data into:
        *   `selected_candidates.csv` (Top N for Oracle)
        *   `unlabeled_pool.csv` (The rest for inference)
5.  **Oracle**:
    *   Runs FEP (Free Energy Perturbation) on the top candidates to get ground truth labels.
    *   **Or**: Runs MM-GBSA (Approximate Free Energy) for faster, cheaper validation.
6.  **Trainer (The Brain)**:
    *   Trains an XGBoost Proxy Model on Oracle Results + RDKit descriptors.
    *   Predicts scores for the entire `unlabeled_pool`.
7.  **Elite Selection**:
    *   Selects the best candidates from the pool for the next cycle (80% Best Predicted + 20% Random).
8.  **Redocking**:
    *   Converts Elite SMILES $\to$ 3D SDFs (Required for Pocket2Mol).
9.  **Fine-Tuning Launcher**:
    *   Submits a Vertex AI Custom Training Job to fine-tune Pocket2Mol on the redocked SDFs.

## Task Specifications
The pipeline logic relies on **JSON Task Templates** in `workspace/`.
*   These JSONs are "Patched" at runtime by `pipeline_e2e.py` to inject dynamic arguments (like `loop_id`).
*   Key Specs: `*task-spec.json`.
