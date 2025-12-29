# Component: The Oracle (Ground Truth)

The Oracle provides the high-fidelity signal that trains the Proxy Model.

## 1. Production Oracle: GROMACS FEP
Free Energy Perturbation (FEP) is the gold standard for computational affinity prediction.

### Mechanics
*   **Method:** Alchemical Transformation (coupling/decoupling the ligand from the solvent and protein).
*   **Protocol:** Non-Equilibrium Fast Growth or Equilibrium Sampling (BAR).
*   **Complexity:** The system runs **13 Lambda Windows** per molecule to ensure convergence.

### Infrastructure (The Heavy Lift)
*   **Compute:** Google Cloud Batch.
*   **Scale:** Each molecule spawns 13 parallel GPU tasks.
*   **Cost:** ~$83.00 per molecule (Standard L4 GPU).
*   **Optimization:** Future roadmap includes **Spot Instance** support via checkpoint-resume to reduce cost to ~$25.

### Data Flow
1.  **Input:** Ligand SDF + Protein PDB.
2.  **Setup:** `fep_setup.py` generates topology (`topol.top`) and structure (`system.gro`).
3.  **Run:** 13 parallel simulations (Minimization $\to$ Equilibration $\to$ Production).
4.  **Analyze:** Bennett Acceptance Ratio (BAR) calculates $\Delta G$.

---

## 2. Development Oracle: Mock FEP
To allow rapid iteration of the pipeline logic without bankruptcy, a synthetic oracle is used.

*   **Logic:** $Score = \text{Docking} + \text{Noise} + \text{QED\_Penalty}$.
*   **Role:** Validates the software architecture (BigQuery, Training, Loop) for $0 cost.
*   **Switch:** Controlled by the pipeline definition (`al_pipeline_def.py`).
