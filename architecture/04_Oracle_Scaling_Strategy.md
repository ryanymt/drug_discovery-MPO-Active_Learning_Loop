# Architecture: Scaling Strategy (MM-GBSA & Batch)

**Context**: Production scoring of 100,000 molecules using MM-GBSA.
**Scale**: 1000 - 100,000 Tasks.
**Workload**: ~1 hour per task (Molecular Dynamics simulation).

## 1. The Constraints

Scaling to 100,000 concurrent simulations presents specific infrastructure challenges driven by **Spot Availability** and **API Rate Limits**.

### A. Quota & Hardware Availability
*   **Resource**: NVIDIA L4 GPUs.
*   **Constraint**: Obtaining thousands of L4 GPUs in a single zone simultaneously is difficult, particularly for Spot instances. But this can be mitigated by using **Multi-Region Dispatch**, or **multi-specification** dispatch.
*   **Solution**: **Multi-Region Dispatch** involves configuring the Batch Job to allow multiple zones (e.g., `us-central1`, `us-east1`, `europe-west4`), enabling the scheduler to hunt for available capacity globally. **Multi-specification** dispatch involves submitting multiple jobs with different specifications (e.g., `n1-highcpu-16`, `g2-standard-12`), enabling the scheduler to hunt for available capacity of different types.

### B. The "Startup Storm"
*   **Scenario**: Launching 10,000 VMs simultaneously creates a spike in API calls and download requests.
*   **Impact**:
    *   **Container Registry**: Throttling on `docker pull`.
    *   **GCS**: 429 (Too Many Requests) errors on input downloads.
*   **Mitigation**: Native Batch Throttling (see below).

---

## 2. Execution Strategy: "Many Small Tasks"

For the previously tested FEP workload during development(7+ hours), a "Checkpoint & Resume" strategy was required to handle preemptions. However, the move to MM-GBSA (~1 hour runtime) allows for a simpler, more robust approach.

### The Strategy
**"Fail Fast, Retry Often"**

1.  **Atomicity**: Each task processes a single molecule (or a small batch of 10).
2.  **No Checkpointing**: The simulation is treated as atomic. If a node is preempted at minute 55, the task is marked `FAILED` and requeued from minute 0.
3.  **Rationale**:
    *   **Simplicity**: Eliminates complex watchdog scripts and synchronization logic.
    *   **Efficiency**: Losing 1 hour of compute on a Spot instance is financially negligible compared to the engineering cost of maintaining stateful resume logic.
    *   **Throughput**: Google Cloud Batch handles the retry loop automatically.

---

## 3. Infrastructure Implementation

### Native Batch Throttling
To mitigate the "Startup Storm" and manage Quota usage without an external orchestrator, the `parallelism` field in the Batch Job Specification is utilized.

#### Configuration Example
```json
{
    "taskGroups": [
        {
            "taskCount": 100000,
            "parallelism": 500
        }
    ]
}
```

#### Behavior
1.  **Queueing**: All 1000 (or 100,000) tasks are submitted to the Google Cloud Batch queue immediately.
2.  **Throttling**: The service ensures only 100 (or 500) VMs are active at any given second.
3.  **Rolling Execution**: As one task completes (or fails), a new task is dequeued to take its place.
4.  **Benefit**: This flattens the API usage curve and keeps the active GPU count within the project's quota limits, preventing "Quota Exceeded" errors.

### Artifact Management
To avoid performance degradation from listing valid outputs (The "Small File Problem"):
*   **No Listing**: Downstream steps do not scan the output directory for `*.csv`.
*   **Manifest Injection**: The job is passed a pre-computed CSV Manifest.
*   **Determinism**: Output paths are deterministic (`shard_${TASK_ID}.csv`).
*   **Consolidation**: A "Map-Reduce" style script aggregates the results by iterating through the known task IDs rather than traversing the filesystem.
