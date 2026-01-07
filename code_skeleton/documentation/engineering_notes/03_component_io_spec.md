# 03. Component I/O Specifications

This document defines the exact Input and Output paths for every component in the pipeline.
**Base Path**: `gs://drug-discovery-mvp-docking-results` (Referred to as `BASE`).
**Dynamic Path**: `.../output/{loop_id}/...` (All components namespace output by `loop_id`).

## 1. Generation (`pocket2mol`)
*   **Input**: None (for Cycle 1) or `checkpoint.pt` (Future Cycles).
    *   *Note*: Requires `models/pretrained.pt` and target PDB (currently `4yhj.pdb`) from `gs://drug-discovery-mvp-input-data`.
*   **Output**: `BASE/output/{loop_id}/generation/SMILES.txt`
    *   Format: Raw text, one SMILES string per line.

## 2. Filtering (Parallel)
### 2a. Gnina
*   **Input**: `BASE/output/{loop_id}/generation/SMILES.txt`
*   **Output**: `BASE/input/{loop_id}/filters/gnina_out` (Directory)
    *   Format: Directory containing SDF files and logs. *Note: Path includes `input` in logic, likely legacy naming. Check `pipeline_e2e.py`.*

### 2b. TxGemma
*   **Input**: `BASE/output/{loop_id}/generation/SMILES.txt`
*   **Output**: `BASE/input/{loop_id}/filters/txgemma_results.csv`
    *   Format: CSV with columns `smiles`, `toxicity`, `affinity`.

### 2c. RDKit
*   **Input**: `BASE/output/{loop_id}/generation/SMILES.txt`
*   **Output**: `BASE/input/{loop_id}/filters/rdkit_scores.csv`
    *   Format: CSV with `smiles`, `qed`, `logp`, `sa`.

## 3. Joiner (`join_results.py`)
*   **Inputs**: All outputs from Step 2.
*   **Output**: `BASE/output/{loop_id}/selection/joined_results.csv`
    *   Format: Consolidated CSV with all scores.

## 4. Selection (`select_active_learning.py`)
*   **Input**: `BASE/output/{loop_id}/selection/joined_results.csv`
*   **Outputs**:
    1.  `BASE/output/{loop_id}/selection/selected_candidates.csv` (Top N for Oracle)
    2.  `BASE/output/{loop_id}/selection/unlabeled_pool.csv` (Remaining candidates)

## 5. Oracle
### Mock FEP (`mock_fep_batch.py`)
*   **Input**: `BASE/output/{loop_id}/selection/selected_candidates.csv`
*   **Output**: `BASE/output/{loop_id}/oracle/fep_results.csv`

### Prod FEP (GROMACS)
*   **Input**: `BASE/output/{loop_id}/selection/selected_candidates.csv`
*   **Output**: `BASE/output/{loop_id}/oracle/fep_results.csv`
    *   Note: Real FEP may produce complex artifacts; result CSV extracts the final `delta_g`.

### Prod MM-GBSA (GROMACS)
*   **Input**: `BASE/output/{loop_id}/selection/selected_candidates.csv`
*   **Output**: `BASE/output/{loop_id}/oracle/mmpbsa_results.csv`
    *   **Logic**: `run_mmpbsa.sh` iterates over candidates, runs 1ns MD + MM-GBSA, and aggregates results.

## 6. Trainer (`train_xgboost.py`)
*   **Inputs**:
    1.  `BASE/output/{loop_id}/oracle/fep_results.csv` (Labels)
    2.  `BASE/output/{loop_id}/selection/unlabeled_pool.csv` (Features)
*   **Outputs**:
    1.  `BASE/output/{loop_id}/training/proxy_model.json` (Artifact)
    2.  `BASE/output/{loop_id}/training/proxy_predictions.csv` (Scores for Pool)

## 7. Elite Selection (`select_elite.py`)
*   **Input**: `BASE/output/{loop_id}/training/proxy_predictions.csv`
*   **Output**: `BASE/output/{loop_id}/selection_v2/elite_candidates.csv`

## 8. Redocking (`redock_batch.py`)
*   **Input**: `BASE/output/{loop_id}/selection_v2/elite_candidates.csv`
*   **Output**: `BASE/output/{loop_id}/training_data/sdfs` (Directory of SDFs)

## 9. Fine-Tuning (`finetune_launcher.py`)
*   **Input**: `BASE/output/{loop_id}/training_data/sdfs`
*   **Output**: `BASE/output/{loop_id}/models/finetuned` (Model Artifacts)
