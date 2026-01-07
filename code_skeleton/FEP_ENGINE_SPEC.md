# FEP Production Engine Specification

## 1. Objective
Create a Python module (`fep_runner.py`) to execute the **production phase** of the Free Energy Perturbation (FEP) workflow. This module picks up where `fep_setup.py` leaves off.

**Input:** Equilibrated system files (`.gro`, `.top`, `.ndx`).
**Output:** A single floating-point value: `delta_G` (kcal/mol), plus logs.

---

## 2. Technical Requirements

### A. Lambda Windowing
The script must run a series of simulations at different "alchemical" states (lambda values).
*   **Standard Scheme:** 20 windows minimum for smooth convergence.
    *   Coulomb (Electrostatics): `0.0, 0.1, ..., 1.0`
    *   VdW (Van der Waals): `0.0, 0.05, ..., 1.0` (often done sequentially or simultaneously depending on soft-core potential settings).
*   **Parallelism:**
    *   *Option A (Sequential):* One GPU runs Window 0 -> Window 1 -> ... (Too slow).
    *   *Option B (Parallel):* The script should utilize GROMACS' internal MPI or simply loop if running on a single powerful node.
    *   *Cloud Batch Context:* Since we are provisioning 1 L4 GPU per "Ligand Job", this script will likely run the windows **sequentially** on that single GPU, OR (better) the Cloud Batch job should have been sharded by *Window* instead of by *Ligand*.
    *   *Constraint:* For this specification, assume **Sequential execution on 1 GPU** for simplicity, or basic thread-parallelism.

### B. Simulation Steps (Per Window)
For each $\lambda$ window:
1.  **Energy Minimization:** `gmx grompp` -> `gmx mdrun` (steepest descent).
2.  **NVT Equilibration:** 100ps, restraining heavy atoms.
3.  **NPT Equilibration:** 100ps, pressure coupling.
4.  **Production MD:** 5ns (or user configurable), saving `dhdl.xvg` energy differences.

### C. Analysis (BAR Method)
After all windows complete:
1.  Collect all `dhdl.xvg` files.
2.  Run `gmx bar` (Bennett Acceptance Ratio).
3.  Parse the output to extract:
    *   `Delta G` (DG)
    *   `Standard Error` (Err)

---

## 3. Implementation Interface (Python)

The engineers should implement a class `FEPSimulator` with the following contract:

```python
import subprocess
from typing import List, Tuple, Dict

class FEPSimulator:
    def __init__(self, work_dir: str, topology_file: str, coordinate_file: str):
        self.work_dir = work_dir
        self.top = topology_file
        self.gro = coordinate_file
        # Standard lambda schedule (example)
        self.lambdas = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

    def _prepare_mdp(self, lambda_val: float, step: str) -> str:
        """
        Generates the .mdp configuration file for a specific window/step.
        Must enable 'free-energy = yes' and set 'init-lambda-state'.
        """
        pass

    def run_window(self, window_idx: int, lambda_val: float):
        """
        Executes EM -> NVT -> NPT -> PROD for a single window.
        """
        # 1. Generate MDPs
        # 2. Run grompp & mdrun (using subprocess)
        # 3. Handle GPU binding (-nb gpu)
        pass

    def run_all_windows(self):
        """
        Iterates through self.lambdas and calls run_window().
        """
        for idx, lam in enumerate(self.lambdas):
            print(f"Starting Window {idx} (Lambda={lam})...")
            self.run_window(idx, lam)

    def analyze_results(self) -> Tuple[float, float]:
        """
        Runs 'gmx bar' on the resulting .xvg files.
        Returns (delta_G_kcal_mol, error_estimate).
        """
        cmd = ["gmx", "bar", "-f", "production_*.xvg", ...]
        # ... parse logic ...
        return delta_g, error

if __name__ == "__main__":
    # Integration with the Pipeline
    sim = FEPSimulator(work_dir="./work", topology_file="topol.top", coordinate_file="eq.gro")
    sim.run_all_windows()
    dg, err = sim.analyze_results()
    print(f"FEP_RESULT: DeltaG={dg}, Error={err}")
```

## 4. Key GROMACS Flags for Cloud
*   `-nb gpu`: Ensure non-bonded interactions run on the L4.
*   `-pme gpu`: If possible, offload PME to GPU.
*   `-update gpu`: Offload integration (requires modern GROMACS).
*   `-ntomp <N>`: Set OpenMP threads to match the vCPU count (e.g., 4).

## 5. Deliverables
1.  `fep_runner.py`: The Python module described above.
2.  `mdp_templates/`: Folder containing template `.mdp` files (`em.mdp`, `nvt.mdp`, `prod.mdp`) with placeholders for lambda values.
