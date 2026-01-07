# AI-Driven Drug Discovery MVP: Target-to-Hit Workflow on Google Cloud

## Project Overview

This project provides a Minimum Viable Product (MVP) for an automated, scalable, cloud-native workflow for the early stages of drug discovery. The primary goal is to computationally screen potential drug candidates against a biological target by modeling the target's natural flexibility.

The workflow leverages containerization with Docker and high-performance computing on Google Cloud, specifically using Google Cloud Batch and Vertex AI Pipelines to execute and orchestrate computationally intensive scientific tasks like molecular dynamics and virtual screening.

## The Vision: A Scalable, Active Learning Drug Discovery Platform

This project aims to demonstrate a next-generation "ChemOps" platform built on Google Cloud, showcasing four key capabilities:

1.  **Iterative Active Learning Loop:**
    *   **Goal:** Move beyond static screening to dynamic optimization.
    *   **Method:** We use `Pocket2Mol` to generate novel candidates, screen them with a multi-parameter filter (RDKit, Gnina, TxGemma), and then use a "Ground Truth" physics-based feedback loop (GROMACS FEP) to train a proxy model (XGBoost). This loop continuously improves the quality of generated molecules over time.

2.  **Vertex AI MLOps & Lineage:**
    *   **Goal:** Automate complex scientific workflows with full traceability.
    *   **Method:** By orchestrating the entire lifecycle (Generation -> Filtering -> Simulation -> Training) on **Vertex AI Pipelines**, we ensure reproducibility, automated metadata tracking (Lineage), and seamless model management.

3.  **Massive Scalability with Google Cloud Batch:**
    *   **Goal:** Drastically reduce time-to-insight for computationally intensive simulations.
    *   **Method:** We leverage **Google Cloud Batch** to burst to thousands of cores (CPUs & GPUs) on demand. This is demonstrated by running highly parallelized Free Energy Perturbation (FEP) calculations (e.g., 13+ concurrent GPU nodes per candidate) to replace weeks of sequential lab work with hours of cloud compute.

4.  **Unlimited Analytics with BigQuery:**
    *   **Goal:** Handle millions of compounds and simulation results without bottlenecks.
    *   **Method:** All pipeline results—from docking scores to ADMET properties and FEP energies—are streamed directly into **BigQuery**. This allows for real-time querying, comprehensive analytics, and the ability to scale data storage infinitely without managing database infrastructure.

---

## The Completed MVP Workflow

The end-to-end process takes a single protein structure and identifies potential "hit" compounds through a multi-stage computational funnel, orchestrated as a cohesive pipeline.

### 1. Target Preparation & Ensemble Generation

This step moves beyond a static "lock-and-key" model by simulating the protein's natural movement and flexibility.

*   **Process:** A Molecular Dynamics (MD) simulation is run to generate thousands of different conformations of the target protein over time.
*   **Tool:** `GROMACS`
*   **Outcome:** A trajectory file containing a "movie" of the protein's movement.

### 2. Conformation Clustering & PDBQT Conversion

The trajectory file is analyzed to select a diverse set of structures that represent the protein's most common shapes.

*   **Process:** A Python script clusters the conformations, and the representative structure from each cluster is converted to the PDBQT format required for docking.
*   **Tools:** `Python`, `MDAnalysis`, `Open Babel`
*   **Outcome:** A small ensemble of representative protein structures in `.pdbqt` format.

### 3. High-Throughput Virtual Screening (Docking)

Each representative protein structure is used to screen a library of potential drug molecules (ligands).

*   **Process:** Each ligand is computationally "docked" into the binding site of each protein structure, and the software calculates a binding affinity score.
*   **Tool:** `GNINA` (a fork of AutoDock Vina)
*   **Outcome:** A ranked list of ligands based on their predicted binding scores.

---

## Getting Started

### Prerequisites

*   A Google Cloud project with the following APIs enabled:
    *   Batch API
    *   Vertex AI API
    *   Artifact Registry API
    *   Cloud Storage API
*   A service account with the necessary IAM permissions.
*   Docker installed locally.
*   Google Cloud SDK installed and configured.

### Setup

> [!IMPORTANT]
> Please refer to [CONTAINER_GUIDELINES.md](CONTAINER_GUIDELINES.md) for important rules regarding Dockerfile best practices (e.g., NO ENTRYPOINT).

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/your-username/your-repo-name.git
    cd your-repo-name
    ```
2.  **Authenticate with Google Cloud:**
    ```bash
    gcloud auth login
    gcloud config set project YOUR_PROJECT_ID
    ```
3.  **Build and push the Docker containers:**
    *   **GROMACS:**
        ```bash
        gcloud builds submit --region YOUR_REGION --config gromacs-container/cloudbuild.yaml
        ```
    *   **Dataproc:**
        ```bash
        gcloud builds submit --region YOUR_REGION --config dataproc-container/cloudbuild.yaml
        ```
    *   **GNINA:**
        ```bash
        gcloud builds submit --region YOUR_REGION --config gnina-container/cloudbuild.yaml
        ```
4.  **Create the Cloud Storage buckets:**
    ```bash
    gsutil mb gs://drug-discovery-mvp-input-data
    gsutil mb gs://drug-discovery-mvp-md-trajectories
    gsutil mb gs://drug-discovery-mvp-docking-results
    ```

### Cloud Storage Buckets Explained

The workflow relies on three Cloud Storage buckets to manage data between steps. Here's what they are and what they should contain:

*   **`gs://drug-discovery-mvp-input-data`**: This bucket holds all the initial input files required to start the pipeline.
    *   **Contents**:
        *   `1aki.pdb`: The initial protein structure.
        *   GROMACS `.mdp` files (`ions.mdp`, `minim.mdp`, `nvt.mdp`, `production.mdp`): Configuration files for the molecular dynamics simulation.
        *   `cluster.py`: The Python script for the clustering and conversion step.
        *   `config.txt`: The configuration file for the GNINA docking step.

*   **`gs://drug-discovery-mvp-md-trajectories`**: This bucket stores the output from the GROMACS molecular dynamics simulation.
    *   **Contents**:
        *   `production.xtc` and `production.gro`: The trajectory and final structure files from the simulation. These are used as input for the clustering step.

*   **`gs://drug-discovery-mvp-docking-results`**: This bucket stores the intermediate and final results of the workflow.
    *   **Contents**:
        *   `representative_structures/`: This directory is created by the `dataproc` step and contains the clustered PDBQT files.
        *   `gnina_output/`: This directory is created by the `gnina` step and contains the final docking scores and poses.
        *   This bucket also serves as the `PIPELINE_ROOT` for Vertex AI, storing pipeline artifacts.

### Running the Workflow

The entire workflow is orchestrated using Vertex AI Pipelines. To run the pipeline, you will need to:

1.  **Define the pipeline in Python.** (A sample script can be found in the `notebooks` directory).
2.  **Compile the pipeline to a JSON file.**
3.  **Submit the pipeline to Vertex AI.**

A Jupyter Notebook, `demo_notebook.ipynb`, has been created to demonstrate and visualize the results from each stage of the completed pipeline. This notebook allows for an interactive exploration of:

*   The molecular dynamics trajectory.
*   The clustered, representative protein structures.
*   The final docking scores and poses of the screened ligands.
