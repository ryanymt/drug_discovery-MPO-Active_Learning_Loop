# Component: Pocket2Mol Generator

## 1. Overview
Pocket2Mol is an E(3)-Equivariant Graph Neural Network (GNN) designed for structure-based drug design. It generates 3D ligands atom-by-atom inside a protein binding pocket.

*   **Role:** The "Actor" (Generates candidates).
*   **Input:** Protein Pocket Structure (`.pdb`).
*   **Output:** 3D Ligand Structure (`.sdf`).

## 2. Operational Modes

### A. Inference (Sampling)
*   **Command:** `sample_for_pdb.py`
*   **Logic:**
    1.  Loads the protein pocket.
    2.  Autoregressively places atoms (C, N, O...) and bonds based on the learned probability distribution.
    3.  Refines geometry.
*   **Compute:** Runs on **Cloud Batch** (CPU or GPU). Can scale to 100,000 molecules easily.

### B. Fine-Tuning (Retraining)
*   **Command:** `train.py` (Patched for Resume).
*   **Logic:** "Weighted Retraining" strategy.
    1.  Loads a pre-trained checkpoint (`pretrained.pt`).
    2.  Loads a new "Elite Dataset" (Top 5% of previous cycle).
    3.  Updates weights to maximize likelihood of generating the Elite molecules.
*   **Compute:** Runs on **Vertex AI Custom Training** (Requires GPU, e.g., Nvidia T4/L4).

## 3. Data Integration
*   **Input Data:** Requires a directory of `.sdf` files aligned to the pocket.
*   **Pre-Processing:** The `PocketLigandPairDataset` class converts SDFs into an LMDB cache for fast training.
*   **Output Model:** Saves checkpoints (`1000.pt`) to GCS.

## 4. Architectural Notes
*   **Equivariance:** The model understands rotation/translation. Rotating the protein results in the generated ligand rotating accordingly.
*   **Chemistry Compliance:** While powerful, it can generate invalid valencies (e.g., Nitrogen with 5 bonds). The system relies on **RDKit Sanitization** downstream to filter these out.
