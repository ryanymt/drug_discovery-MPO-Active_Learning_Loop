# Target Plan: A Two-Phase Approach to AI-Driven Drug Discovery

## Current Status and Next Steps

**Current Status:**
*   The individual components of a potential drug discovery pipeline (GROMACS, Dataproc, GNINA, Pocket2Mol, Txgemma) have been analyzed.
*   A critical bug in the GNINA component's output script has been identified and fixed.
*   Initial documentation outlining the flow and target architecture has been created.
*   The **TxGemma container is not functional**, which blocks the creation of the complete Phase 1 pipeline. The model files are also currently ignored by version control, preventing a reproducible build.

**Immediate Next Steps:**

1.  **Fix the TxGemma Container:** This is the highest priority.
    *   **Action:** Modify the `.gitignore` file to stop ignoring the `txgemma-container/container/model/` directory. This is essential for building a self-contained, reproducible container image.
    *   **Action:** Debug the `txgemma-container/container/Dockerfile` and `predict.py` script to ensure the container builds successfully and the prediction script runs without errors. This includes installing all necessary dependencies and verifying the model can be loaded correctly.

2.  **Develop the BigQuery Loader Component:**
    *   **Action:** Create a new containerized component responsible for reading CSV results from GCS (from GNINA and TxGemma) and loading them into a predefined BigQuery table.

3.  **Begin Scaling Pocket2Mol and GNINA:**
    *   **Action:** Start drafting the necessary modifications for the `pocket2mol` and `gnina` batch job specifications to handle large-scale, parallel execution as required for Phase 1.

---

## Phase 1: The "Linear" Foundation

**Objective:** Create a high-throughput, scalable pipeline for generating and screening a large number of molecules, focusing on engineering performance and reproducibility.

**Required Work:**

1.  **Component Creation: ADMET Scoring (TxGemma)**
    *   **Task:** Develop a new pipeline component for ADMET prediction using what you've referred to as "TxGemma".
    *   **Implementation:**
        *   Create a new container image with the necessary libraries to run the TxGemma model.
        *   Write a script that takes a batch of molecules (e.g., in SDF or SMILES format) as input, runs inference, and outputs scores.
        *   Define a new task specification JSON (e.g., `txgemma-task-spec.json`) in the `workspace/` directory for this component.

2.  **Pipeline Orchestration: New Vertex AI Pipeline**
    *   **Task:** Create a new pipeline definition file (e.g., `workspace/phase1_pipeline.py`) that orchestrates the linear workflow.
    *   **Workflow:**
        1.  **Pocket2Mol:** Generate ~100,000 molecules.
        2.  **GNINA:** Dock all generated molecules against the target protein.
        3.  **TxGemma:** Run ADMET scoring on all generated molecules.
        4.  **BigQuery Load:** Collect results from GNINA and TxGemma and load them into a structured BigQuery table.

3.  **Component Scaling & Modification**
    *   **Pocket2Mol:** The current `pocket2mol-batch-job.json` is for a single run. This needs to be adapted to generate a massive batch and handle the output efficiently, likely sharding the results into multiple files on GCS.
    *   **GNINA:** The existing gnina-task is designed for a few representative structures. This step must be heavily parallelized to handle ~100,000 docking calculations. We will likely use a Batch job with a large number of parallel tasks, where each task docks a subset of the molecules.

4.  **Data Infrastructure: BigQuery Setup**
    *   **Task:** Design and create a BigQuery table to store the results.
    *   **Schema:** The table should include columns for molecule ID, SMILES string, docking score (GNINA), ADMET properties (from TxGemma), and any other relevant metadata.
    *   **Implementation:** A new component (or a script within an existing one) will be needed to format the data from GCS and load it into BigQuery.

---

## Phase 2: The "Iterative" Active Learning Loop

**Objective:** Transform the linear pipeline into a closed-loop system that uses high-fidelity simulation feedback to iteratively improve the molecule generation model.

**Required Work:**

1.  **New Component Development: Acquisition & Training**
    *   **Acquisition Selector:** A new component that queries the BigQuery results from a Phase 1 run and applies a selection strategy (e.g., top-k, uncertainty sampling) to choose the best ~100 molecules for the next step.
    *   **Pocket2Mol Trainer:** A new component responsible for fine-tuning the Pocket2Mol model. It will take the 100 molecules and their corresponding high-fidelity scores (from GROMACS FEP) as input and output a new, updated model checkpoint to GCS.

2.  **High-Fidelity Validation: GROMACS (FEP)**
    *   **Task:** Enhance the existing GROMACS component or create a new one to perform Free Energy Perturbation (FEP) calculations.
    *   **Implementation:** This is significantly more complex than the current MD simulation. It will involve a more detailed setup, much longer run times, and require powerful GPU resources. The output will be a precise binding energy score for each of the 100 selected molecules.

3.  **Pipeline Architecture for Active Learning**
    *   **Task:** Re-architect the Vertex AI pipeline to support iteration and state management.
    *   **Implementation:** The pipeline will be designed as a loop.
        *   It will start with an initial Pocket2Mol model.
        *   Each iteration will execute the Phase 1 Engine, the Acquisition Selector, the GROMACS FEP validation, and the Pocket2Mol Trainer.
        *   The key challenge is managing the model state: the trainer component in loop N must produce a model checkpoint that is then used by the generation component in loop N+1.

4.  **Metrics and Tracking**
    *   **Task:** Implement robust tracking to generate the "Hero Plot".
    *   **Implementation:** The BigQuery table will be extended to track which iteration a molecule was generated in. A separate process or notebook will query this data to plot the improvement of the Multi-Parameter Optimization (MPO) score over time, proving the value of the active learning loop.