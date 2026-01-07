# Component: Scoring Engines (The Filter Layer)

To optimize molecules, we must measure them. We use a multi-tiered scoring strategy.

## 1. Gnina (Docking Affinity)
*   **Role:** The primary "Affinity" signal.
*   **Technology:** Deep Learning-based Molecular Docking (CNN). Based on AutoDock Vina.
*   **Metric:**
    *   **Vina Affinity:** Predicted binding energy (kcal/mol). Lower is better (e.g., -9.0).
    *   **CNN Score:** Probability that the pose is correct (0-1).
*   **Infrastructure:** Runs on **Cloud Batch** (GPU accelerated).
*   **Input:** Ligand `.sdf`, Receptor `.pdb`.

## 2. RDKit (Chemoinformatics)
*   **Role:** The "sanity check" and property calculator.
*   **Technology:** Rule-based chemical library.
*   **Metrics:**
    *   **QED:** Quantitative Estimate of Drug-likeness (0-1). Measures solubility, size, etc.
    *   **SA:** Synthetic Accessibility (1-10). How hard is it to make?
    *   **LogP:** Lipophilicity.
*   **Infrastructure:** Runs locally or on lightweight Batch nodes (CPU).

## 3. TxGemma (Safety/ADMET)
*   **Role:** The "Safety" filter.
*   **Technology:** Large Language Model (LLM) fine-tuned for toxicology.
*   **Metrics:** Toxicity probability, LD50 prediction.
*   **Infrastructure:** Requires GPU inference. Currently implemented but calibration is ongoing.

## 4. Parallel Scoring Architecture
Instead of scoring sequentially, we launch massive parallel jobs on **Google Cloud Batch**.
*   **Batch Size:** 1,000 - 100,000 molecules.
*   **Scripting:** We use unified Bash scripts (`run_gnina_loop2.sh`) to iterate through files on GCS mounts, avoiding data transfer bottlenecks.
