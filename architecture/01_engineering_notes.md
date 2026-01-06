# 99. Engineering Log: Production System Verification (Phase 1)

**Status**: Active Learning Loop - Cycle 1 Complete / Cycle 2 In Progress
**Scale**: Cycle-1 100,000 Molecules, Cycle-2 10,000 Molecules
**Date**: January 2026

## Section 1. Executive Summary
This document records the engineering details, configuration, and observations from the "Phase 1" production scale-up. The objective was to validate the full end-to-end active learning loop, generating 100,000 molecules, scoring them with high-fidelity physics, training the surrogate models for the cycle-2. In cycle-2, generate 10,000 molecules using the re-trained generator, and screen them with the same pipeline to compare the scores.

A high-affinity binder (**-66.25 kcal/mol**) was successfully identified, and a Proxy Oracle was trained with **R²=0.82**.


### Compute Cost
**Compute cost of running one round** (from baseline generation to pocket2mol fine-tuning): **~850 USD** 

*   **Spot Instances**: Used for 100% of the stateless Batch execution to reduce costs by ~70%.
*   **Fault Tolerance**: A "Many Small Tasks" strategy is employed. If a Spot node is preempted, only a small shard (e.g., 1,000 molecules) is retried, minimizing wasted compute. Checkpoint mechanism can be implemented to get better cost efficiency.
*   **Containerization**: All compute steps are encapsulated in Docker containers, pushed to Artifact Registry

### Performance & Scalability Analysis
**Total Runtime**: ~15 Hours (Cycle 1 Generation to Cycle 2 Initiation)
*   **Current Hardware**: Limited to **100 GPUs** (L4/T4 mix).
*   **Bottleneck**: The Oracle (MM-GBSA) step takes ~10 hours for 1000 molecules with this concurrency limit.
*   **Projection**: With **500 GPUs**, the entire loop runtime is estimated to drop to **~3 Hours**, enabling rapid daily iteration cycles.

### Tasks breakdown are as below. 

#### Cloud Batch Jobs

| Component | CPU | Memory | GPU | Tasks | Parallelism | Duration | Notes |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Pocket2Mol** | 4 | 16GB | L4 | 100 | 50 | 1h 48 min | Generate baseline 100k molecules, ~70k molecules yield |
| **RDkit** | 4 | 16GB | None | 11 | 1 | 6 min | RDKit scoring on 100k molecules |
| **Gnina** | 16 | 16GB | T4 | 100 | 100 | 28 min | Docking score on 100k molecules. CPU heavy tasks, almost zero GPU usage |
| **TxGemma** | 4 | 16 | L4 | 100 | 100 | 11 mins | Toxicity score on 100k molecules |
| **Gromacs** | 16 | 48GB | L4 | 1000 | 100 | 9h 51 min | MM-GBSA score on 1000 molecules using 1ns MD. CPU heavy, only 40% GPU load  |
| **Proxy Model Prediction (XGboost)** | 8 | 32GB | None | 1 | 1 | 3 min | Predcit Delta G score on 99,000 molecules using proxy model |

#### Vertex AI Training
| Component | CPU | Memory | GPU | Duration | Notes |
| :--- | :--- | :--- | :--- | :--- | :--- | 
| **Pocket2Mol Fine-Tuning** | 12 | 85GB | A100 40GB | 1 | 1 | 2h 45 min | Pocket2Mol fine-tuning on 10k molecules |


---
## Section 2 : Baseline (Cycle 1)

## Step I: Large Scale Generation (Pocket2Mol)
A distributed batch job was executed to generate 100,000 molecules using the `Pocket2Mol` generator. To optimize for Spot Instance reliability, a "Many Small Tasks" strategy was employed.

*   **Job ID**: `prod-pocket2mol-100k-v2`
*   **Scale**: 100 tasks × 1,000 molecules (50 GPUs concurrent)
*   **Output**: 70,506 Valid Molecules (70% Validity Rate)

### Key Assets
*   **Spec**: `workspace/prod-pocket2mol-100k-v2.json`
*   **Artifacts**: SMILES text files and 3D SDF structures stored in `gs://ryanymt/output/generated/shard_*/`.

To manage the massive file count, a **Manifest Consolidation** step (`consolidate_manifest.py`) was implemented, aggregating all 70k+ paths into a single CSV for downstream ingestion. 

### BigQuery
**bioops_platform** dataset was created at this point to store the molecules and all the scoring values generated from following tasks will be updated here. **BigQuery Schema** is detailed at the end of this document.


---

## Step II: The Screening Funnel
The 100k candidates were passed through a multi-stage screening filter to remove toxic, insoluble, or non-drug-like compounds before expensive physics scoring.

### A. RDKit Filtering
QED (Quantitative Estimation of Drug-likeness), SA Score (Synthetic Accessibility), and LogP were calculated. This step ran efficiently on standard CPUs.

*   **Job**: `prod-rdkit-100k` (10 CPUs)
*   **Throughput**: 100k molecules processed in < 15 minutes.

### B. Toxicity Screening (TxGemma)
The `google/txgemma-9b-predict` model was utilized to flag potential clinical toxicity.

*   **Job**: `prod-txgemma-100k` (100 x L4 GPUs)
*   **Outcome**: Flagged molecules with `toxicity_label=1` were excluded from the Oracle candidates.

### C. Molecular Docking (Gnina) - *Challenges Encountered*
The CNN-based docking engine `gnina` was deployed to estimate binding affinity.

> [!TIP]
> **Performance Optimization**
> Gnina was found to be CPU-bound (Vina search). A standard GPU node leaves the GPU 95% idle.
> **Fix**: Use `n1-highcpu-16` (16 vCPUs) and run 16 concurrent workers per T4 GPU. This saturates the hardware and increases throughput from 4 to **30 molecules/minute**.

---

## Step III: The Oracle (MM-GBSA)
For the "Ground Truth" confirmation, top candidates were selected to undergo Molecular Mechanics with Generalized Born Surface Area (MM-GBSA) calculations using GROMACS.

### Selection Strategy
**876 Candidates** were selected for this expensive step:
1.  **Exploitation (Top 80% )**: Best docking scores.
2.  **Exploration (Random 20% )**: Diversity sampling.
3.  **Safety**: Must be Non-Toxic (`TxGemma=0`) and Drug-Like (`QED>0.4`).

### Execution Details
*   **Job**: `prod-mmpbsa-100k`
*   **Compute**: 100 x `g2-standard-12` (L4 GPU).
*   **Throughput**: ~1 hour per molecule (1.0 ns Simulation).

### Results
*   **Success Rate**: 92% (810/876 completed).
*   **Top Hit**: **-66.25 kcal/mol** (Significant affinity).
*   **Population Mean**: -35.5 kcal/mol.

---

## Step IV: Closing the Loop (Learning)
With the Ground Truth data established, the surrogate models were trained to guide the next generation cycle.

### Surrogate Model Training (XGBoost)
A "Proxy Oracle" was trained to predict MM-GBSA DeltaG scores based solely on molecular fingerprints (ECFP4), bypassing the need for docking poses.

*   **Training Set**: 810 labeled molecules.
*   **Validation**:
    *   **R² Score**: **0.82** (Excellent predictivity).
    *   **RMSE**: 5.20 kcal/mol.
*   **Inference**: The model scored the remaining 70,000 molecules in seconds (`prod-inference-100k`).


### Elite Selection 
To prevent "Mode Collapse" (generating similar molecules repeatedly), a clustering-based selection strategy was implemented for the next training set.

*   **Strategy**:
    1.  Filter for valid properties ( `SA < 4.0` and `QED > 0.5` (Only drug-like molecules)).
    2.  Cluster the top 20k candidates (Morgan Fingerprints) based on affanity.
    3.  Select Centroids (Tier 1: highest score in each structure cluster) + Random Diversity (Tier 2- random backfilled to reach 10k).
*   **Outcome**: Created a dataset of **10,000 Elite Candidates** (`cycle-02-baseline`).

---

### Step V: Retraining of Generator (Cycle 2 Prep)
The loop continues with the retraining of the generator.

**Fine-Tuning**: Job `cycle-02-pocket2mol-finetune-a100-v5` successfully retrained Pocket2Mol on the Elite Dataset.
    *   **Config**: `max_iters: 25000` (10 Epochs), `batch_size: 4`, `use_apex: False`.
    *   **Vertex AI Metadata**: Logs execution, input artifacts (SDF files), and output artifacts (Model PT).

---

## Section 3: Cycle 2 

**Cycle 2 Generation**: Job `prod-pocket2mol-cycle2-10k` runs to generate 10,000 new molecules from the biased distribution using re-trained pocket2mol-v2.
**Scoring**: RDKit, TxGemma, Proxy Model Prediction (XGboost) are used to score the newly generated 10,000 molecules

**Results** are published in **result** directory.


---
## Section 4: BigQuery Schema Reference 

### 1. `molecule_registry` (Source of Truth)
Contains the unique identity of every generated molecule across all cycles.
*   **`molecule_hash`** (STRING, PK): SHA256 of SMILES.
*   **`smiles`** (STRING): Canonical SMILES string.
*   **`global_id`** (STRING): Original generation ID (e.g., `gen_shard0_idx1`).
*   **`sdf_gcs_path`** (STRING): Path to 3D structure (e.g., `gs://.../0.sdf`).
*   **`run_id`** (STRING): Batch identifier (e.g., `cycle1_100k`, `cycle2_10k`).
*   **`created_at`** (TIMESTAMP): Insertion time.

### 2. `screening_results` (Fast Scores)
Contains property predictions and docking scores.
*   **`molecule_hash`** (STRING, FK): Links to Registry.
*   **`docking_score`** (FLOAT): Gnina CNN Affinity.
*   **`toxicity_label`** (INTEGER): 0 (Safe) or 1 (Toxic) from TxGemma.
*   **`qed_score`** (FLOAT): Quantitative Estimation of Drug-likeness (0-1).
*   **`sa_score`** (FLOAT): Synthetic Accessibility Score (1-10).
*   **`final_score`** (FLOAT): Unified affinity score (Ground Truth or Proxy).
*   **`score_method`** (STRING): Source of score ('gromacs', 'xgboost').

### 3. `final_affinity` (Ground Truth)
Contains the expensive, high-fidelity MM-GBSA calculations.
*   **`molecule_hash`** (STRING, FK): Links to Registry.
*   **`final_deltag`** (FLOAT): Binding Free Energy (kcal/mol).

### 4. `elite_candidates_v1` (Active Learning Pool)
Snapshot of top candidates selected for retraining.
*   **`molecule_hash`** (STRING, FK): Links to Registry.
*   **`smiles`** (STRING): For easy access.
*   **`final_score`** (FLOAT): Score used for ranking.
*   **`selection_strategy`** (STRING): 'centroid' or 'diversity'.
*   **`dataset_type`** (STRING): 'cycle_2_baseline'.
