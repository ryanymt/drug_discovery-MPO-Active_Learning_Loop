# Pocket2Mol Fine-Tuning Strategy

## Context
In our Active Learning loop (`01-prod-pipeline.json`), we improve the model quality by iteratively fine-tuning the "Generator" (Pocket2Mol) on the highest-scoring candidates discovered by the Oracle (MM-GBSA).

## Approaches Considered

### Option A: Reinforcement Learning (RL)
*   **Verdict**: Rejected.
*   **Reason**: High instability and cost (requires online Oracle interaction).

### Option B: Conditional Generation (Control Codes)
*   **Verdict**: Rejected.
*   **Reason**: Requires massive property-labeled datasets that we don't have.

### Option C: Iterative Fine-Tuning (Selected Strategy)
Retrain the model weights on a small, curated dataset of "Elite" candidates ("Hill Climbing").
*   **Upside**: Simple, robust, offline.
*   **Risks**: Catastrophic Forgetting, Mode Collapse.

## Production Implementation (Cycle 2 Status)

We successfully implemented **Option C** with specific optimizations for scale ("Weighted Retraining").

### 1. Data Preparation: The "Retrieval" Strategy
Unlike the initial plan to "Redock" SMILES (which causes coordinate drift), we implemented a **Retrieval Strategy**.
*   **Logic**: Since Pocket2Mol is a 3D generative model, it generated the original specific 3D conformer that yielded the high score.
*   **Mechanism**: We trace the `molecule_hash` back to the original `shard_ID` and `index` in the Generation output (`gs://ryanymt/output/generated/`) and retrieve the *exact* original SDF.
*   **Start**: `select_elite.py` (BigQuery -> CSV).
*   **Bridge**: `prep_p2m_data.py` (CSV -> Tarball of SDFs).

### 2. Training Infrastructure: A100 Acceleration
We migrated from P100s to **NVIDIA A100 GPUs** for the retraining phase.
*   **Why**: The fine-tuning dataset (10k molecules) is small, but the graph processing (PyG) is data-intensive.
*   **Optimization**: We implemented **Multiprocessing LMDB Generation** to speed up data loading from 60 mins -> 4 mins.
*   **Config**: `batch_size: 4`, `max_iters: 25000` (10 Epochs).

### 3. Diversity Controls
To mitigate Mode Collapse, we employ a "Robust Recipe" for the training set:
*   **80% Centroids**: Representative high-scorers from structural clusters.
*   **20% Exploration**: Random diverse samples.
*   **Result**: The fine-tuned model achieved a **99.9% Hit Rate** in Cycle 2, proving the strategy works.
