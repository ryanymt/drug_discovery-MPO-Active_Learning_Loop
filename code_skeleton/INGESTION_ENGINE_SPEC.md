# Engineering Directive: BigQuery Data Ingestion Implementation

## 1. Objective
Refactor the `join-container` to transition from a file-based aggregator (CSV merging) to a **BigQuery-based Ingestion Engine**. The goal is to separate data by "Fidelity" into the three tiers defined in the architecture: GOLD (Physics), SILVER (Docking), and BRONZE (ML Predictions).

## 2. Architecture Context
*   **Current State:** `join.py` reads individual result files from GCS, merges them into one CSV, and uploads it back to GCS.
*   **Target State:** `join.py` should read result files and **stream** them into specific BigQuery tables.
*   **Why?** This enables the "Active Learning Loop" to query historical data instantly without parsing thousands of CSVs.

## 3. Implementation Requirements

### A. Dependency Update
*   Update `join-container/Dockerfile` or `requirements.txt` to include:
    *   `google-cloud-bigquery`
    *   `google-cloud-storage` (already present)
    *   `pandas` (for easier data manipulation)

### B. Logic Refactoring (`join.py`)
The script must accept a new argument `--data_tier` [GOLD, SILVER, BRONZE] and `--bq_dataset`.

#### 1. Ingestion Strategy (The "Router")
Instead of a single logic path, implement three distinct handlers:

**Case 1: Tier = GOLD (FEP Results)**
*   **Input:** `.../mol_*/analyze/output/results.json` (from `run_fep_workflow.py` --analyze).
*   **Schema Mapping:**
    *   `molecule_hash` -> SHA256(SMILES)
    *   `delta_g` -> `float`
    *   `uncertainty` -> `float`
    *   `trajectory_uri` -> `gs://.../mol_*/run/` (Construct this path)
*   **Target Table:** `bioops_platform.validation_results`

**Case 2: Tier = SILVER (Gnina Results)**
*   **Input:** `.../gnina_output/*.csv` (or JSON).
*   **Schema Mapping:**
    *   `molecule_hash` -> SHA256(SMILES)
    *   `cnn_score` -> `float`
    *   `affinity` -> `float`
    *   `pose_uri` -> `gs://.../docked_poses/`
*   **Target Table:** `bioops_platform.screening_results`

**Case 3: Tier = BRONZE (Proxy Predictions)**
*   **Input:** `.../predictions/*.csv` (from XGBoost).
*   **Schema Mapping:**
    *   `predicted_delta_g` -> `float`
    *   `model_version` -> (Input Argument)
*   **Target Table:** `bioops_platform.proxy_predictions`

#### 2. The UPSERT Logic (Crucial)
*   **Problem:** We might run FEP on the same molecule twice (e.g., in Loop 1 and Loop 5).
*   **Requirement:** We want the *latest* high-fidelity result.
*   **Approach:** Use BigQuery `MERGE` statement instead of simple `INSERT`.
    ```sql
    MERGE bioops_platform.validation_results T
    USING Staging S
    ON T.molecule_hash = S.molecule_hash
    WHEN MATCHED THEN
      UPDATE SET delta_g = S.delta_g, run_id = S.run_id
    WHEN NOT MATCHED THEN
      INSERT (...) VALUES (...)
    ```
    *Implementation Note:* The script can load data into a temporary table first, then run the MERGE query.

### C. Interface Changes
Update the `argparse` definition to support:
```python
parser.add_argument("--data_tier", choices=["GOLD", "SILVER", "BRONZE"], required=True)
parser.add_argument("--bq_dataset", default="bioops_platform")
parser.add_argument("--run_id", help="Vertex Pipeline Run ID for lineage")
```

## 4. Verification
The engineer should verify success by running a test job and querying:
```sql
SELECT count(*) FROM `lifescience-project-469915.bioops_platform.validation_results` WHERE run_id = 'test-run-001'
```
