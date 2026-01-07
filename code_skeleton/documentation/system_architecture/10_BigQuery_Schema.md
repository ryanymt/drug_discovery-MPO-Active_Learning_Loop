# BigQuery Schema Definition

## 1. Dataset: `bioops_platform`

## 2. Table: `molecule_registry`
**Description:** The master index of all generated molecules.
**Partitioning:** By `run_id`.

```json
[
  {"name": "molecule_hash", "type": "STRING", "mode": "REQUIRED", "description": "SHA256 of Canonical SMILES (Primary Key)"},
  {"name": "smiles", "type": "STRING", "mode": "REQUIRED", "description": "Canonical SMILES"},
  {"name": "global_id", "type": "STRING", "mode": "NULLABLE", "description": "Unique Generated ID (e.g. gen_shard0_idx1)"},
  {"name": "sdf_gcs_path", "type": "STRING", "mode": "NULLABLE", "description": "GS Path to 3D Structure (e.g. gs://.../0.sdf)"},
  {"name": "run_id", "type": "STRING", "mode": "REQUIRED", "description": "Production Batch ID (e.g. prod_100k_v1)"},
  {"name": "created_at", "type": "TIMESTAMP", "mode": "NULLABLE"}
]
```

## 3. Table: `screening_results`
**Description:** Phase 2 Scores (Docking, RDKit, AI).
**Clustering:** By `molecule_hash`.

```json
[
  {"name": "molecule_hash", "type": "STRING", "mode": "REQUIRED", "description": "Foreign Key to registry"},
  {"name": "docking_score", "type": "FLOAT", "mode": "NULLABLE", "description": "Gnina CNN Score / Vina Affinity (Lower is better usually, depending on metric)"},
  {"name": "qed_score", "type": "FLOAT", "mode": "NULLABLE", "description": "Quantitative Estimation of Drug-likeness (0-1)"},
  {"name": "sa_score", "type": "FLOAT", "mode": "NULLABLE", "description": "Synthetic Accessibility Score (1-10)"},
  {"name": "logp", "type": "FLOAT", "mode": "NULLABLE"},
  {"name": "mw", "type": "FLOAT", "mode": "NULLABLE", "description": "Molecular Weight"},
  {"name": "toxicity_label", "type": "INTEGER", "mode": "NULLABLE", "description": "TxGemma Prediction (0=Safe, 1=Toxic)"}
]
```

## 4. Table: `final_affinity`
**Description:** The "Golden" binding score, sourced from either expensive Physics (Oracle) or ML Proxy.
**Granularity:** One row per molecule per run.

```json
[
  {"name": "molecule_hash", "type": "STRING", "mode": "REQUIRED"},
  {"name": "affinity_score", "type": "FLOAT", "mode": "NULLABLE", "description": "Delta G (Negative is better). Source depends on 'method'."},
  {"name": "method", "type": "STRING", "mode": "REQUIRED", "description": "Enum: 'MMGBSA', 'XGBOOST', 'FEP'"},
  {"name": "uncertainty", "type": "FLOAT", "mode": "NULLABLE", "description": "Standard Deviation (for Physics) or Variance (for ML)"},
  {"name": "stage", "type": "STRING", "mode": "NULLABLE", "description": "e.g. 'oracle_v1', 'proxy_v1'"}
]
```

## 5. SQL View: `v_leaderboard`
**Logic:** A unified view to rank molecules by their best available score.

```sql
SELECT
    r.smiles,
    a.affinity_score,
    a.method,
    s.docking_score,
    s.qed_score,
    a.molecule_hash
FROM `bioops_platform.molecule_registry` r
JOIN `bioops_platform.final_affinity` a ON r.molecule_hash = a.molecule_hash
LEFT JOIN `bioops_platform.screening_results` s ON r.molecule_hash = s.molecule_hash
ORDER BY a.affinity_score ASC
```

## 4. Table: `oracle_results` (MM-GBSA / FEP)
**Description:** Expensive physics validation for the top 1%.

```json
[
  {"name": "molecule_hash", "type": "STRING", "mode": "REQUIRED"},
  {"name": "method", "type": "STRING", "mode": "REQUIRED", "description": "MMPBSA or FEP"},
  {"name": "delta_g", "type": "FLOAT", "mode": "NULLABLE"},
  {"name": "uncertainty", "type": "FLOAT", "mode": "NULLABLE"},
  {"name": "status", "type": "STRING", "mode": "NULLABLE", "description": "SUCCESS or FAILURE"}
]
```

## 5. Table: `proxy_predictions`
**Description:** ML predictions for the molecules we didn't run physics on.

```json
[
  {"name": "molecule_hash", "type": "STRING", "mode": "REQUIRED"},
  {"name": "predicted_delta_g", "type": "FLOAT", "mode": "NULLABLE"},
  {"name": "model_version", "type": "STRING", "mode": "NULLABLE"}
]
```

## 6. SQL View: `v_comprehensive_leaderboard`
**Logic:** Combines Real Physics (Gold) with AI Predictions (Bronze) to rank everything.

```sql
CREATE OR REPLACE VIEW `bioops_platform.v_comprehensive_leaderboard` AS
SELECT
    r.smiles,
    r.molecule_hash,
    COALESCE(o.delta_g, p.predicted_delta_g) AS best_estimate_score,
    CASE
        WHEN o.delta_g IS NOT NULL THEN 'GOLD_PHYSICS'
        ELSE 'BRONZE_AI'
    END AS confidence_tier,
    s.qed_score,
    s.toxicity_score
FROM `bioops_platform.molecule_registry` r
LEFT JOIN `bioops_platform.oracle_results` o ON r.molecule_hash = o.molecule_hash
LEFT JOIN `bioops_platform.proxy_predictions` p ON r.molecule_hash = p.molecule_hash
LEFT JOIN `bioops_platform.screening_results` s ON r.molecule_hash = s.molecule_hash
```
