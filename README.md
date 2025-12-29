# Cloud-Native Active Learning for Drug Discovery
### Author
* **Author:** Ryan Ye Min Thein
* **GitHub:** [ryanymt](https://github.com/ryanymt)

![Pareto Front](results/plots/pareto_front_final.png)

## Overview
This repository documents the architecture and results of a **Closed-Loop AI Drug Discovery Platform** built on Google Cloud. 

By orchestrating Generative AI (Pocket2Mol), Physics-Based Simulation (FEP), and Surrogate Modeling (XGBoost), the platform achieved a **12.8x enrichment** in hit rates compared to random generation.

### The Active Learning Concept
![Generic Active Learning Loop](diagrams/Active_Learning_loop_generic.png)

## Key Results
*   **Hit Rate:** Increased from 2.7% (Baseline) to **35.0%** (Cycle 2).
*   **Affinity:** Shifted population mean by **+0.7 kcal/mol**.
*   **Scale:** Validated architecture for 100,000+ molecules using Google Cloud Batch.

## Architecture Highlights
The system implements a specific MPO (Multi-Parameter Optimization) strategy using the following components:

![ChemOps Platform Architecture](diagrams/Drug_Discovery-MPO-Active_Learning.png)

The system is built on a "Lakehouse" pattern using:
*   **Vertex AI Pipelines:** For end-to-end orchestration.
*   **Google Cloud Batch:** For massive parallel docking and FEP simulations (160+ GPUs).
*   **BigQuery:** For structured chemical data analytics.

## Documentation
Dive into the detailed system design:

1.  [Architecture Overview](architecture/00_Architecture_Overview.md) - High-level diagrams and flow.
2.  [Generative Engine](architecture/01_Pocket2Mol_Generator.md) - Pocket2Mol fine-tuning process.
3.  [The Oracle](architecture/03_Oracle_FEP.md) - Using GROMACS FEP as Ground Truth.
4.  [The Critic](architecture/04_Proxy_Model_XGBoost.md) - Surrogate modeling for rapid screening.
5.  [Strategy: Diversity](architecture/07_Strategy_Mode_Collapse.md) - Mitigating mode collapse in generative models.

## Results
See the full breakdown of the Cycle 1 vs Cycle 2 analysis in the [Final Project Report](results/final_project_report.md).

---
*Note: This repository contains architectural documentation and result artifacts. Development code and container definitions are omitted.*

