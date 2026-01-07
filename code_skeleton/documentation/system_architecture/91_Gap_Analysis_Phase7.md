# Gap Analysis: Closing the Active Learning Loop

## Current State
The current pipeline (`00-mock-pipeline.json` / `01-prod-pipeline.json`) is **Open Loop**:
`Generation` $\to$ `Screening` $\to$ `Selection` $\to$ `Oracle (FEP)` $\to$ **[STOP]**

## Goal State (Closed Loop)
`Generation` $\to$ ... $\to$ `Oracle` $\to$ **`Training (XGBoost)`** $\to$ **`Scoring (Proxy)`** $\to$ **`Selection (Elite + Random + Replay)`** $\to$ **`Fine-Tuning (Pocket2Mol)`** $\to$ **`Loop 2 Generation`**

## Identified Gaps (The "Missing Half")

### Gap 1: The "Unlabeled Pool" Handover
*   **Problem**: Selector picks the Top N for Oracle and discards the rest.
*   **Gap**: Need to save the remaining molecules as `unlabeled_pool.csv` for the Proxy Model to score.

### Gap 2: XGBoost Training Integration
*   **Problem**: Training was manual.
*   **Target**: Pipeline Component taking `oracle_results` + `screening_results` $\to$ `Model Artifact`.

### Gap 3: The "Elite Selection" Logic
*   **Problem**: Manual selection.
*   **Target**: Component taking `unlabeled_pool` + `model` $\to$ Inference $\to$ Diversity Filter $\to$ `elite_candidates.csv`.

### Gap 4: The "SMILES-to-3D" Bridge (Redocking)
*   **Problem**: Fine-tuning requires 3D SDFs; pipeline passes SMILES.
*   **Target**: Component taking `elite_candidates.csv` $\to$ Gnina/Obabel $\to$ `training_data.tar.gz`.

### Gap 5: Fine-Tuning Trigger
*   **Problem**: Long-running job risks timeout in standard container.
*   **Target**: Component that submits a **Vertex AI Custom Training Job** and waits/polls for completion.

### Gap 6: Cycle 2 Generation
*   **Problem**: Manual launch.
*   **Target**: Parameterize Generation Component to accept optional `checkpoint.pt`.

## Strategic Roadmap

### Phase 1: Data Preservation (The Bridge)
*   [ ] Modify Selector to output:
    1.  `candidates_for_oracle.csv` (Top N)
    2.  `candidates_for_inference.csv` (The Rest)

### Phase 2: The "Brain" Components
*   [ ] Build Trainer Component (XGBoost).
*   [ ] Build Inference Component (Batch Predict).
*   [ ] Build Elite Selector (Logic + Diversity).

### Phase 3: The "Actor" Update
*   [ ] Build Redocker Component (SMILES $\to$ SDF).
*   [ ] Build Fine-Tuning Launcher (Vertex Custom Job).
*   [ ] Wire Checkpoint into Generation Component.


## Phase 10: 100k Production Scale-Up Analysis (Completed Jan 2026)

We have successfully executed the full "Closed Loop" at scale (100k Generation).

### 1. Scale Validation
*   **Generation**: Scaled Pocket2Mol to 100k molecules (50 GPUs).
*   **Oracle**: Scaled MM-GBSA to 876 molecules (100 L4 GPUs).
*   **Throughput**: Achieved ~20 mins/molecule for Oracle, proving L4 viability over FEP.

### 2. Active Learning Validation
We proved that the loop actually works (i.e., the model learns).

*   **Cycle 1 (Baseline)**:
    *   **Population**: 70,000 Unique Molecules.
    *   **Affinity Mean**: -7.7 kcal/mol.
    *   **Elite Rate (<-20)**: ~23%.
    *   **Top Hit**: -66 kcal/mol (Found 1 "Needle in Haystack").

*   **Cycle 2 (AI-Driven)**:
    *   **Population**: 10,000 Generated Molecules.
    *   **Affinity Mean**: **-36.9 kcal/mol**.
    *   **Improvement**: **4.8x shift** in binding energy distribution.
    *   **Elite Rate (<-20)**: **>99%**.
    *   **Conclusion**: The fine-tuned model (trained on just ~800 data points) successfully learned to generate high-affinity candidates exclusively.

### 3. Remaining Gaps (Optimizations)
1.  **Gnina Coordinates**: The "Redocking" gap (Gap 4) was partially solved by retrieving original SDFs, but Gnina drift remains a risk. Future: Use .
2.  **Cost**: MM-GBSA is cheap (/bin/bash.50), but scaling to 1M molecules requires further reduction (e.g., shorter MD, implicit solvent only).
3.  **Diversity**: Cycle 2 is heavily focused (Mode Collapse). Future: Increase "Exploration" weight (Random Sampling) in elite selection.
