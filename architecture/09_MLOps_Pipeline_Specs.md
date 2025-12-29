# MLOps: Pipeline Architecture & Orchestration

## 1. Overview
The platform uses **Vertex AI Pipelines (Kubeflow)** to orchestrate the multi-stage drug discovery workflow. The platform maintains two distinct versions of the pipeline to balance development speed with production accuracy.

## 2. Pipeline Variants

### A. Development Pipeline (`00-mock-pipeline.json`)
*   **Purpose:** Rapid iteration and logic verification.
*   **Oracle:** Uses the **Synthetic Oracle (Python)** to simulate FEP results.
*   **Cost:** Minimal (Standard CPU/GPU compute for filtering).
*   **Speed:** End-to-end execution in ~30 minutes.

### B. Production Pipeline (`01-prod-pipeline.json`)
*   **Purpose:** Generating high-fidelity leads for lead optimization.
*   **Oracle:** Uses **GROMACS FEP (GCP Batch)** with 13 parallel lambda windows per molecule.
*   **Cost:** High (~$83 per molecule).
*   **Speed:** End-to-end execution in ~12 hours.

---

## 3. Orchestration Mechanics: "Spec Patching"
The current implementation uses a "Patching" pattern to handle the dynamic nature of Batch Job specifications inside a static Kubeflow graph.

1.  **Template Spec:** A hardcoded JSON string representing a Google Cloud Batch job.
2.  **Patch Component:** A Python-based component that uses Regex to inject runtime variables (e.g., `loop_id`, `gcs_path`, `molecule_count`) into the JSON.
3.  **Submit Component:** A generic Batch Job submission component that takes the patched JSON and executes it via the GCP Batch API.

---

## 4. Identified Gaps & Roadmap

### Gap 1: The "Open" Loop
Currently, the pipelines stop at the **Oracle** stage. The data is ingested into BigQuery, but the subsequent **XGBoost Training** and **Pocket2Mol Fine-Tuning** are currently manual or external triggers.

*   **Roadmap:** Integrate the Vertex AI Custom Training Job as a formal component at the end of the DAG to achieve "Full-Loop Automation."

### Gap 2: Spot Instance Risk in Prod
The `01-prod-pipeline.json` uses Spot instances for the 7-hour GROMACS windows without a checkpoint-resume mechanism.

*   **Roadmap:** Implement the "Watchdog" syncer logic in the production container to enable safe use of Spot instances, reducing FEP costs by 70%.

### Gap 3: Data Lineage (Metadata)
While the pipeline runs, manual paths in GCS are used. 

*   **Roadmap:** Fully leverage **Vertex AI Metadata** to pass artifacts (SDFs, CSVs) between components using `Input[Artifact]` and `Output[Artifact]` types instead of hardcoded GCS strings.
