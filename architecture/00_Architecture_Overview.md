# System Architecture: AI-Driven Drug Discovery Platform

## 1. Executive Summary
This platform implements a closed-loop **Active Learning** system for de novo drug design. It combines generative AI (Pocket2Mol) with physics-based validation (FEP) and surrogate modeling (XGBoost) to iteratively explore chemical space and optimize molecules for binding affinity, drug-likeness, and safety.

The system is fully **Cloud-Native**, leveraging Google Cloud Batch for massive parallelism and Vertex AI for machine learning operations.

---

## 2. High-Level Workflow (The Active Learning Loop)

The core logic revolves around a cycle of **Generation**, **Scoring**, and **Learning**.

### Cycle 1: The Baseline (Cold Start)
1.  **Generation:** The pre-trained Pocket2Mol model generates 1,000 diverse molecules binding to the target pocket.
2.  **Scoring:** Molecules are evaluated by RDKit (Properties) and Gnina (Docking).
3.  **Selection:** A diverse subset is chosen for labeling.
4.  **Labeling:** The Oracle (FEP or Mock) assigns high-fidelity scores.
5.  **Proxy Training:** An XGBoost model learns to predict the Oracle score from molecular fingerprints.

### Cycle 2: The Optimization (Fine-Tuned)
6.  **Prediction:** The Proxy Model scores a larger pool of candidates.
7.  **Elite Selection:** The top candidates (High Predicted Affinity + Good Properties) form the "Elite Set".
8.  **Fine-Tuning:** Pocket2Mol is retrained on the Elite Set to shift its distribution.
9.  **Sampling:** The new model generates better molecules.

---

## 3. System Components Diagram

![ChemOps Platform Architecture](../diagrams/Drug_Discovery-MPO-Active_Learning.png)

## 4. Key Infrastructure Decisions

| Component | Technology | Reasoning |
| :--- | :--- | :--- |
| **Orchestration** | Vertex AI Pipelines | Managed Kubeflow, lineage tracking, experiment management. |
| **Compute (HPC)** | Google Cloud Batch | Cost-effective scaling for 100k+ simulations (Gnina/FEP). Handles Spot instances. |
| **Compute (ML)** | Vertex AI Training | Managed environment for PyTorch/XGBoost training jobs. |
| **Storage** | Google Cloud Storage | Blob storage for massive unstructure data (SDF, Trajectories). |
| **Warehouse** | BigQuery | Structured query engine for filtering millions of candidates by score. |

---
###[Notes: Code is not published as part of this repository]
## 5. Directory Map
*   `workspace/`: Operational scripts and temporary artifacts.
*   `pocket2mol-container/`: Source code for the Generator.
*   `gnina-container/`: Source code for the Docker engine.
*   `gromacs-container/`: Source code for the FEP engine.
*   `documentation/`: Architectural records.
