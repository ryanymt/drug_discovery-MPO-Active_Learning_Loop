# Drug Discovery Pipeline Component Documentation

This document outlines the components of the drug discovery pipeline, their inputs, and their outputs.

## GROMACS: Molecular Dynamics Simulation

*   **Purpose**: Performs a molecular dynamics (MD) simulation on a protein to generate an ensemble of conformations, capturing the protein's natural flexibility.
*   **Container**: `gcr.io/lifescience-project-469915/drug-discovery-containers/gromacs-gpu:2023.3.29`
    *   **Dockerfile**: `gromacs-container/Dockerfile`
    *   **Script**: `gromacs-container/run_gromacs.sh`
*   **Batch Job**: `workspace/gromacs-task-spec.json`
*   **Inputs**:
    *   **GCS Bucket**: `gs://drug-discovery-mvp-input-data`
        *   `1aki.pdb`: The initial protein structure.
        *   GROMACS `.mdp` files (`ions.mdp`, `minim.mdp`, `nvt.mdp`, `production.mdp`): Configuration files for the MD simulation.
*   **Outputs**:
    *   **GCS Bucket**: `gs://drug-discovery-mvp-md-trajectories`
        *   `production.xtc`: The trajectory file from the simulation.
        *   `production.gro`: The final structure file from the simulation.
        *   `em.gro`: The energy-minimized structure.

## Dataproc: Clustering and PDBQT Conversion

*   **Purpose**: Clusters the conformations from the MD simulation and converts the representative structures to the PDBQT format required for docking.
*   **Container**: `gcr.io/lifescience-project-469915/drug-discovery-containers/dataproc:1.2.1`
    *   **Dockerfile**: `dataproc-container/Dockerfile`
    *   **Script**: `dataproc-container/cluster.py`
*   **Batch Job**: `workspace/dataproc-task-spec.json`
*   **Inputs**:
    *   **GCS Bucket**: `gs://drug-discovery-mvp-md-trajectories`
        *   `em.gro`: The energy-minimized structure.
        *   `nvt.gro`: The structure after NVT equilibration.
*   **Outputs**:
    *   **GCS Bucket**: `gs://drug-discovery-mvp-docking-results/representative_structures`
        *   `representative_0.pdbqt`, `representative_1.pdbqt`: The representative protein structures in PDBQT format.

## GNINA: Docking

*   **Purpose**: Performs high-throughput virtual screening (docking) of a library of ligands against the representative protein structures.
*   **Container**: `gcr.io/lifescience-project-469915/drug-discovery-containers/gnina:latest`
    *   **Dockerfile**: `gnina-container/container/Dockerfile`
    *   **Script**: `gnina-container/container/run_gnina_docking.sh`
*   **Batch Job**: `workspace/gnina-task-spec.json`
*   **Inputs**:
    *   **GCS Bucket**: `gs://drug-discovery-mvp-docking-results/sharded_smiles_output/sharded_smiles`
        *   `smiles_part_000.txt`: A file containing SMILES strings for the ligands.
    *   **GCS Bucket**: `gs://drug-discovery-mvp-input-data`
        *   `1aq1.pdb`: The protein structure.
*   **Outputs**:
    *   **GCS Bucket**: `gs://drug-discovery-mvp-docking-results/gnina_output`
        *   `predictions_000.csv`: A CSV file containing the docking scores for each SMILES string.
        *   `docking_log_000.txt`: A log file for the docking run.

## Pocket2Mol: Molecule Generation

*   **Purpose**: Generates novel molecules based on a protein's 3D binding pocket. This can be used to create a library of candidate ligands for the GNINA docking step.
*   **Container**: Not explicitly defined in a batch job, but a Dockerfile is provided.
    *   **Dockerfile**: `pocket2mol-container/Dockerfile`
    *   **Scripts**: `pocket2mol-container/sample.py`, `pocket2mol-container/sample_for_pdb.py`
*   **Inputs**:
    *   A protein PDB file.
    *   The center coordinates of the binding pocket.
*   **Outputs**:
    *   Generated molecule files (e.g., SDF).

## Txgemma: Property Prediction

*   **Purpose**: Predicts properties of molecules. This can be used to filter or rank the docked ligands from the GNINA step.
*   **Container**: `gcr.io/lifescience-project-469915/drug-discovery-containers/txgemma:v2`
    *   **Dockerfile**: `txgemma-container/container/Dockerfile`
    *   **Script**: `txgemma-container/container/predict.py`
*   **Batch Job**: `workspace/txgemma-task-spec.json`
*   **Inputs**:
    *   **GCS Bucket**: `gs://drug-discovery-mvp-docking-results/sharded_smiles_output/sharded_smiles`
        *   `smiles_part_000.txt`: A file containing SMILES strings.
*   **Outputs**:
    *   **GCS Bucket**: `gs://drug-discovery-mvp-docking-results/txgemma_output`
        *   `predictions_000.csv`: A CSV file containing the predicted properties for each SMILES string.
