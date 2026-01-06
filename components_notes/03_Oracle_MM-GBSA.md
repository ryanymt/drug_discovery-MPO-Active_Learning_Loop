# Component: The Oracle (Ground Truth)

The Oracle provides the high-fidelity signal that trains the Proxy Model.

## 1. Production Oracle: GROMACS MM-GBSA
The production oracle utilizes **MM-GBSA (Molecular Mechanics-Generalized Born Surface Area)**. This method offers an optimal trade-off between physical accuracy and computational throughput, allowing thousands of candidates to be scored at a fraction of the cost of FEP.

### Mechanics
*   **Method:** End-Point Free Energy Calculation.
*   **Protocol:** Explicit solvent Molecular Dynamics (MD) simulation followed by implicit solvent energy analysis.
*   **Theory:** $\Delta G_{bind} = \Delta H - T\Delta S \approx \Delta E_{MM} + \Delta G_{sol} - T\Delta S$.
*   **Implementation:** The system uses `gmx_MMPBSA` (AmberTools binding) on GROMACS trajectories.

### Infrastructure
*   **Compute:** Google Cloud Batch (Vertex AI Pipelines orchestration).
*   **Scale:** One task per molecule (Embarrassingly Parallel).
*   **Hardware:** NVIDIA L4 GPUs (24GB VRAM).
*   **Cost:** ~$0.20 - $0.50 per molecule (vs $80+ for FEP).
*   **Throughput:** ~15-20 minutes per molecule (1ns simulation).

### Data Flow
1.  **Input:** Ligand SDF (Correct 3D Coordinates) + Protein PDB (AlphaFold/Crystal).
2.  **Topology Generation:** `fep_setup.py` (ACPYPE/AmberTools) generates GROMACS topologies (`.top`, `.gro`).
3.  **Simulation (MD):**
    *   **Minimization:** Steepest descent to remove clashes.
    *   **Equilibration:** NVT/NPT (100ps) to stabilize temperature/pressure.
    *   **Production:** 1ns - 5ns trajectory collection.
4.  **Analysis:** `gmx_MMPBSA` processes the trajectory frames to calculate average Binding Free Energy ($\Delta G$).

---

## 2. Development Oracle: Mock Score
To enable rapid iteration of the pipeline logic without incurring cloud costs, a synthetic oracle is used during development.

*   **Logic:** $Score = \text{Docking} + \text{Noise} + \text{QED\_Penalty}$.
*   **Role:** Validates the software architecture (BigQuery, Training, Loop) for $0 cost.
*   **Switch:** Controlled by the pipeline definition (`al_pipeline_def.py`).
