# 93. MM-GBSA Implementation Strategy

## 1. Goal
Replace the computationally expensive FEP Oracle (13 lambda windows x 1 ns) with a cheaper **MM-GBSA** Oracle (1 standard MD simulation + Post-processing).
*   **Cost Reduction**: ~13x faster simulation time.
*   **Tool**: `gmx_MMPBSA` (Industry standard for GROMACS).

## 2. Container Strategy: `gromacs-mmpbsa`
We will clone `gromacs-container` and modify it.

### Dockerfile Changes
*   **Base**: Keep `gromacs-2023.3` build.
*   **Additions**:
    *   Install `gmx_MMPBSA` via pip (requires `ambertools` which is already in `micromamba`).
    *   Ensure `mpi4py` or necessary deps are present.

```dockerfile
# ... (After micromamba install) ...
RUN micromamba run -n fep_env pip install gmx_MMPBSA
```

## 3. Simulation Strategy (`run_mmpbsa.sh`)
Refactor `run_fep.sh` to remove the lambda loop.

### Workflow
1.  **Setup (`do_setup`)**:
    *   Keep existing `fep_setup.py` (It generates `system.gro`, `topol.top`).
    *   *Note*: FEP setup prepares a dual-topology or specific hybrid structure. For MM-GBSA, we just need a standard complex. The existing setup might be overkill but should still work (Ligand is present).
    *   *Optimization*: Eventually, we can simplify `fep_setup.py` to just `gmx pdb2gmx` + `gmx solvate`, but for now, reusing `fep_setup.py` ensures consistency with the previous pipeline.

2.  **Run (`do_run`)**:
    *   **Minimization**: `minim.mdp`
    *   **Equilibration**: `nvt.mdp` -> `npt.mdp` (Optional, or just NVT)
    *   **Production**: `production.mdp` (1 ns, **Standard MD**, No Free Energy flags).
        *   Output: `production.xtc` (Trajectory), `production.tpr`.

3.  **Analysis (`do_analyze` -> `do_mmpbsa`)**:
    *   **Command**:
        ```bash
        gmx_MMPBSA -O -i mmpbsa.in -cs production.tpr -ct production.xtc -ci index.ndx -cg 1 13
        ```
    *   **Inputs**:
        *   `mmpbsa.in`: Configuration file (see below).
        *   `index.ndx`: Need to generate this to define "Protein" and "Ligand" groups.

## 4. Input Files

### `production_mmgbsa.mdp`
*   Copy `production_fep.mdp`.
*   **Remove**: All `free-energy`, `couple-moltype`, `sc-*` flags.
*   **Keep**: `integrator`, `nsteps`, `pbc`, `tcoupl`.

### `mmpbsa.in`
Create a simple config file.
```text
&general
  sys_name="Prot-Lig-Complex",
  startframe=1,
  endframe=9999,
  interval=10,  # Skip frames for speed
/
&gb
  igb=5,        # GB-OBC(II) model (standard)
  saltcon=0.15,
/
```

## 5. Pipeline Integration
*   **Task Spec**: Create `mmpbsa-task-spec.json`.
*   **Script**: Create `mmpbsa_batch.py` (wrapper around `run_mmpbsa.sh`).
*   **Pipeline**: Add a toggle in `pipeline_e2e.py` to choose `mode="mmpbsa"`.
