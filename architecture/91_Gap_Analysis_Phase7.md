# Gap Analysis: Closing the Active Learning Loop
## Notes: These pipeline gap are implemented of 29 Dec 2025. Keeping this document as development trace

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

