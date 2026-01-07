# 04. Troubleshooting & Gotchas

## Common Issues

### 1. `PipelineParameter` Serialization Error (Type Error)
*   **Symptom**: `TypeError: Object of type PipelineParameterChannel is not JSON serializable` during compilation.
*   **Cause**: Trying to `json.dumps()` a dictionary that contains KFP pipeline parameters (like `loop_id` placeholders) inside the pipeline definition.
*   **Fix**: We use a custom component `patch_spec` (`pipeline_e2e.py`) to inject these parameters at *runtime* inside the container, rather than at compile time. **Do not remove `patch_spec`**.

### 2. Cloud Batch Quota Exceeded
*   **Symptom**: Jobs stuck in `QUEUED` or `SCHEDULED` state for hours.
*   **Check**: Go to IAM & Admin > Quotas. Check `NVIDIA_L4_GPUS` and `NVIDIA_A100_GPUS`.
*   **Fix**: Request quota increase or switch `00-mock-pipeline.json` (CPU-only) for dev work.

### 3. "File Not Found" in Containers
*   **Symptom**: `python: can't open file '.../app/score.py': [Errno 2] No such file or directory`.
*   **Cause**:
    *   Script is missing from GCS `data/scripts/`.
    *   Container entrypoint is overriding the command.
*   **Fix**:
    *   Upload script: `gsutil cp workspace/myscript.py gs://drug-discovery-mvp-docking-results/data/scripts/`.
    *   Verify the `cpgs://... .` command in the Task Spec includes the correct path.

### 4. Pocket2Mol Model Not Found
*   **Symptom**: Generation step fails silently or halts.
*   **Cause**: Image `pocket2mol:latest` expects a pretrained model at `/app/data/pretrained.pt`.
*   **Fix**: Ensure the Docker image was built with the model artifact included, or mount it from GCS.

### 5. GROMACS FEP Failures
*   **Symptom**: `grompp` failures (topology mismatch).
*   **Cause**: Generated molecule topology (from LigParGen/Sobtop) doesn't match the forcefield.
*   **Fix**: Check `01-prod-fep-task-spec.json` inputs. Ensure the `ligand.itp` is generated correctly in the preparation step (hidden inside the Oracle container logic).

## Monitoring
*   **Vertex AI Pipelines**: High-level graph view. Good for seeing which step failed.
*   **Cloud Batch**: Low-level logs. Go here to see `stdout`/`stderr` from the actual python scripts.
*   **Cloud Logging**: Search for `resource.type="batch.googleapis.com/Job"`.
