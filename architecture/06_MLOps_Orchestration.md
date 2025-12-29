# MLOps & Orchestration

The platform is built on Google Cloud's managed services to ensure scalability and reproducibility.

## 1. Vertex AI Pipelines
*   **Role:** The Workflow Engine.
*   **Definition:** Python-based DSL (`al_pipeline_def.py`) compiled to KFP JSON.
*   **Components:** Each step (Generation, Scoring, Training) is a modular containerized component.

## 2. Google Cloud Batch
*   **Role:** The Compute Engine.
*   **Why:** Traditional Kubernetes (GKE) is too expensive/complex for bursty, massive simulations. Batch allows spinning up 1,000 VMs for 1 hour and then shut them down to $0.
*   **Integration:** Python driver scripts submit Batch jobs dynamically using the Python SDK.

## 3. Container Strategy
A suite of specialized Docker images is maintained in Artifact Registry:

*   `pocket2mol-retrain`: Contains PyTorch, PyG, RDKit, and the Generative Model logic.
*   `gnina`: Contains the Gnina binary and OpenBabel.
*   `gromacs-fep`: Contains GROMACS 2023, Acpype, and setup scripts.
*   `proxy-model`: Contains XGBoost, Scikit-Learn, and Google Cloud libraries.

## 4. CI/CD & Script Injection
*   **Pattern:** "Infrastructure as Code, Scripts as Config".
*   **Logic:** The strategy avoids baking rapid-change logic into Docker containers. Instead, driver scripts (`run_gnina_loop2.sh`) are stored in GCS and mounted into the container at runtime.
*   **Benefit:** Allows iterating on logic in seconds without 30-minute container rebuilds.
