# Data Platform Architecture

We employ a "Lakehouse" architecture to manage the diverse data types in drug discovery.

## 1. Storage: Google Cloud Storage (The Lake)
Used for unstructured, heavy, or file-based artifacts.

*   `gs://ryanymt/input`: Immutable references (Protein PDBs, Configs, Scripts).
*   `gs://ryanymt/output`: Pipeline outputs.
    *   `/generated`: Raw shard output from Pocket2Mol (`SMILES.txt`, `.sdf`).
    *   `/prod_100k`: Production run outputs (Gnina logs, TxGemma results).
    *   `/models`: Trained checkpoints (`xgboost_proxy.joblib`, `pocket2mol-v2.pt`).
    *   `/oracle_batches`: Aggregated SDFs for MM-GBSA simulations.

## 2. Warehouse: BigQuery (The House)
Used for structured metrics, querying, and analytics.

**Dataset:** `gcda-apac-sc.bioops_platform`

### Tables (The "Three-Tier" Logic)

#### 1. `molecule_registry` (Source of Truth)
*   **Role**: Identity Management.
*   **Key**: `molecule_hash` (Primary), `run_id` (Partition).
*   **Content**: Mapping of `smiles` to `global_id` and `sdf_gcs_path`.
*   **Logic**: Managed by `upsert_registry.py` (MERGE) to ensure 0 duplicates.

#### 2. `screening_results` (Fast Scores)
*   **Role**: The "Silver" Layer.
*   **Key**: `molecule_hash`.
*   **Content**: `docking_score` (Gnina), `qed_score`, `sa_score`, `toxicity_label`.
*   **Logic**: Ingested via `consolidate_scores.py`.

#### 3. `final_affinity` (Ground Truth)
*   **Role**: The "Gold" Layer (Oracle).
*   **Key**: `molecule_hash`.
*   **Content**: `final_deltag` (MM-GBSA Binding Free Energy).
*   **Logic**: Populated only after expensive validation runs.

#### 4. `elite_candidates_v1` (Active Learning)
*   **Role**: The "Platinum" Layer.
*   **Key**: `molecule_hash`.
*   **Content**: Top 10k candidates selected for Retraining Cycle 2.

## 3. Data Movement Strategy

### The "Upsert" Pattern
To handle the massive redundancy in generative models (where the same molecule is generated thousands of times), we treat BigQuery as an Idempotent Store:
1.  **Generate**: Millions of molecules in GCS shards.
2.  **Consolidate**: Join shards locally or via Dataflow.
3.  **Merge**: Use `upsert_registry.py` to insert only NEW unique hashes into `molecule_registry`, ignoring duplicates.

### The "Unified Score" View
We use SQL logic to present a single `final_score` to downstream consumers:
```sql
COALESCE(gold.final_deltag, silver.docking_score, bronze.proxy_prediction)
```
This allows the Leaderboard to transparently show the "Best Available Truth" for every molecule.
