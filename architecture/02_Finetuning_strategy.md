# Pocket2Mol Fine-Tuning Strategy

## Context
In the Active Learning loop (`01-prod-pipeline.json`), the aim is to improve the quality of generated molecules over successive cycles. The "Generator" (Pocket2Mol) must learn from the high-scoring candidates identified by the Oracle (FEP).

## Approaches Considered

### Option A: Reinforcement Learning (RL)
Directly optimize the generator using the Oracle scores as a reward function (e.g., REINVENT, PPO).
*   **Upside**: Theoretically optimal for maximizing a specific scalar reward (Binding Affinity).
*   **Downside**:
    *   **Sparse Rewards**: High-scoring molecules are rare events.
    *   **Instability**: RL on 3D graph generation is notoriously unstable and hard to tune.
    *   **Complexity**: Requires rewriting the Pocket2Mol loss function and training loop significantly.
    *   **Cost**: Requires online interaction with the Oracle during training (extremely expensive with FEP).

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
    *   **Mode Collapse**: The model may converge to generating only a few specific scaffolds found in the elite set.

## Comparison Table

| Feature | RL | Conditional | Fine-Tuning (Ours) |
| :--- | :--- | :--- | :--- |
| **Implementation Effort** | High | High | **Low** |
| **Stability** | Low | Medium | **High** |
| **Oracle Cost** | High (Online) | Low (Offline) | **Low (Offline)** |
| **Risk** | Non-convergence | Poor Generalization | Mode Collapse |

## Justification for Selection
For the MLOps implementation (Phase 7), **Option C (Iterative Fine-Tuning)** was selected for the following reasons:

1.  **MVP Pragmatism**: The existing `run_finetune_v7.py` script can be reused without modifying the core Pocket2Mol research code.
2.  **Infrastructure Fit**: It maps cleanly to a Batch Pipeline workflow:
    *   `Step N`: Generate Data -> `Step N+1`: Train on Data.
3.  **Risk Management**: "Mode Collapse" is manageable via the **Diversity Filter** in the Selection Step (80% Exploitation / 20% Exploration). Randomness is explicitly injected to prevent the dataset from becoming too narrow.

## Implementation Details (The "SMILES-to-3D" Bridge)
A specific challenge for Pocket2Mol is that it is a **3D Generative Model**—it requires 3D SDF structures as training labels, not just SMILES strings.

**The Pipeline Solution:**
1.  **Selection**: Pick Top N SMILES (`select_elite.py`).
2.  **Redocking**: Use Gnina/RDKit to Generate 3D Poses (`redock_batch.py`).
    *   *Why*: The *predicted* 3D pose of a high-affinity SMILES is treated as the "Ground Truth 3D Structure" for the generator to learn.
3.  **Training**: Run Fine-Tuning on these redocked SDFs (`finetune_launcher.py`).

This "Self-Training" loop allows the 2D-to-3D conversion to happen automatically within the pipeline.
