# Pocket2Mol Fine-Tuning Strategy

## Context
In the Active Learning loop, the aim is to improve the quality of generated molecules over successive cycles. The "Generator" (Pocket2Mol) must learn from the high-scoring candidates identified by the Oracle (FEP / MM-GBSA).

## Approaches Considered

### Option A: Reinforcement Learning (RL)
Directly optimize the generator using the Oracle scores as a reward function (e.g., REINVENT, PPO).
*   **Upside**: Theoretically optimal for maximizing a specific scalar reward (Binding Affinity).
*   **Downside**:
    *   **Sparse Rewards**: High-scoring molecules are rare events.
    *   **Instability**: RL on 3D graph generation is notoriously unstable and hard to tune.
    *   **Complexity**: Requires rewriting the Pocket2Mol loss function and training loop significantly.
    *   **Cost**: Requires online interaction with the Oracle during training (extremely expensive with FEP / MM-GBSA).

### Option B: Conditional Generation (Control Codes)
Train a single model conditioned on property tokens (e.g., `<LogP:High>`, `<Binding:-10>`).
*   **Upside**: One model serves all needs; controllable generation at inference time.
*   **Downside**:
    *   **Data Scarcity**: Requires a massive dataset spanning the *entire* distribution of properties to learn robust embeddings. The high-affinity region is undersampled.
    *   **Architecture**: Requires modifying the Pocket2Mol architecture to accept conditional vectors.

### Option C: Iterative Fine-Tuning (Selected Strategy)
Retrain the model weights on a small, curated dataset of "Elite" candidates ("Hill Climbing").
*   **Upside**:
    *   **Simplicity**: Uses the standard training loop; no architecture changes needed.
    *   **Robustness**: Standard Supervised Learning is stable.
    *   **Offline**: Does not require Oracle interaction during training.
*   **Downside**:
    *   **Catastrophic Forgetting**: The model may "forget" general chemistry rules (valency, ring strain) if trained too long on a small dataset.
    *   **Mode Collapse**: The model may converge to generating only a few specific scaffolds found in the elite set. Experience replay has to be implemented to minimise this. 

## Comparison Table

| Feature | RL | Conditional | Fine-Tuning (**) |
| :--- | :--- | :--- | :--- |
| **Implementation Effort** | High | High | **Low** |
| **Stability** | Low | Medium | **High** |
| **Oracle Cost** | High (Online) | Low (Offline) | **Low (Offline)** |
| **Risk** | Non-convergence | Poor Generalization | Mode Collapse |

## Justification for Selection
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
