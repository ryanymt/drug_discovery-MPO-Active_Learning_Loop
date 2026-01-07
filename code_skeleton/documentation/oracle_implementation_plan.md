# Oracle Implementation Plan: Mock vs. Real

To balance cost, speed, and fidelity during the development of the Active Learning Loop, we utilize a dual-strategy for the "Oracle" component (the source of ground truth).

## 1. Development Mode: The Synthetic Oracle (Mock FEP)

**Objective:** Validate the software architecture (Pipeline, BigQuery, Training, Fine-Tuning) effectively for **$0 cost**.

### The Mechanism
Instead of running GROMACS, we use a lightweight Python script (`mock_fep_batch.py`) that acts as a simulator.

*   **Logic:**
    $$ \text{Score}_{\text{FEP}} = \text{Score}_{\text{Gnina}} + \mathcal{N}(0, \sigma) + \text{Penalty}_{\text{QED}} $$
*   **Why it works:**
    *   **Correlation:** It maintains the physical intuition that "Good Docking usually implies Good FEP".
    *   **Variance:** The noise term $\mathcal{N}$ simulates the fact that Docking is imprecise, forcing the Proxy Model to learn a mapping that isn't just a copy of Gnina.
    *   **Speed:** Labels 1,000 molecules in seconds.

### Implementation Details
*   **Script:** `workspace/mock_fep_batch.py`
*   **Input:** `selected_candidates.csv` (Must contain `docking_score`).
*   **Output:** Generates `fep_final.csv` in the exact GCS path structure required by the pipeline (`fep_batch/mol_{i}/output/`).
*   **Usage:**
    ```bash
    python3 workspace/mock_fep_batch.py --input selected_candidates.csv --count 1000
    ```

---

## 2. Production Mode: High-Fidelity Physics (GROMACS FEP)

**Objective:** Generate chemically accurate Binding Free Energy ($\Delta G$) labels to drive real drug discovery.

### The Mechanism
We use **Google Cloud Batch** to orchestrate massive parallel Molecular Dynamics simulations.

*   **Engine:** GROMACS 2023 (Custom Container).
*   **Method:** Free Energy Perturbation (FEP) / Bennett Acceptance Ratio (BAR).
*   **Scale:** Each molecule requires 13 parallel GPU tasks (Lambda Windows).
*   **Cost:** ~$83.00 per molecule (Standard/On-Demand). (13 tasks x 11 hours x $0.58/hr)

### Implementation Details
*   **Driver:** `workspace/run_fep_batch.py`
*   **Container:** `gromacs-fep:v5` (Includes `run_fep.sh`, `fep_setup.py`).
*   **Infrastructure:**
    *   **Quota:** Requires substantial NVIDIA L4 GPU quota (160+ for parallel batches).
    *   **Data:** Outputs `production.xtc` (trajectory) and `fep_final.csv` (score).
*   **Usage:**
    ```bash
    python3 workspace/run_fep_batch.py --count 10 --id_map candidate_map.csv
    ```

## 3. Transition Strategy

1.  **Phase 1 (Now):** Use **Mock Oracle** to debug the "Outer Loop" (Generation $\to$ Selection $\to$ Training $\to$ Fine-Tuning). Prove the loop creates "Better Mock Scores" over time.
2.  **Phase 2 (Pilot):** Run **Real FEP** on a small batch (10 molecules) to validate the physics engine mechanics (file transfers, formats).
3.  **Phase 3 (Production):** Switch to Real FEP for the full campaign when budget allows.
