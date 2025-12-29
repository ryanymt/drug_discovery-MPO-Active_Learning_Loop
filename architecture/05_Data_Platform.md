# Data Platform Architecture

A "Lakehouse" architecture is employed to manage the diverse data types in drug discovery.

## 1. Storage: Google Cloud Storage (The Lake)
Used for unstructured, heavy, or file-based artifacts.

*   `[INPUT_BUCKET]`: Immutable references (Protein PDBs, Configs).
*   `[RESULTS_BUCKET]`: Pipeline outputs.
    *   `/generated`: SDF files from Pocket2Mol.
    *   `/fep_batch`: Trajectories (`.xtc`), Energy files (`.xvg`).
    *   `/models`: Trained checkpoints (`.pt`, `.json`).

## 2. Warehouse: BigQuery (The House)
Used for structured metrics, querying, and analytics.

**Dataset:** `[PROJECT_ID].chemops_platform`

### Tables
1.  **`validation_results` (Gold Tier)**
    *   **Source:** Oracle (FEP).
    *   **Schema:** `smiles`, `delta_g` (float), `run_id`, `uncertainty`.
    *   **Use:** Training Ground Truth.

2.  **`screening_results` (Silver Tier)**
    *   **Source:** Gnina, RDKit.
    *   **Schema:** `smiles`, `docking_score`, `qed`, `sa`.
    *   **Use:** Selection strategy input.

3.  **`proxy_predictions` (Bronze Tier)**
    *   **Source:** XGBoost.
    *   **Schema:** `smiles`, `predicted_delta_g`.
    *   **Use:** Selecting the "Elite Set" for fine-tuning.

## 3. Data Movement
*   **Joiners:** Custom Python scripts (`join_fep_results.py`) parse raw files from GCS and perform streaming inserts into BigQuery.
*   **Lineage:** All tables include a `run_id` to trace back to the specific pipeline execution.
