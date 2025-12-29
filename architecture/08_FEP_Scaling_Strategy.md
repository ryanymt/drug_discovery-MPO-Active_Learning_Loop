# Architecture: Scalable FEP Strategy (The 13,000 GPU Problem)

**Notes**: The development test was run with 10 molecules out of top 100, not 1,000.
**Context**: Running FEP for 1,000 molecules.
**Scale**: 1,000 Mols $\times$ 13 Windows = **13,000 GPU Tasks**.

## 1. The Constraints

Running 13,000 concurrent GPU tasks is not a code problem; it's a **Physics & Quota** problem.

### A. Quota Wall
*   **Availability**: Getting 13,000 L4 GPUs in a single GCP zone (e.g., `us-central1-a`) is impossible, especially for Spot.
*   **Solution**: **Multi-Region Dispatch**. The Batch Job must allow `allowedLocations: ["us-central1", "us-east1", "us-west1", "europe-west4"]`.

### B. Storage Bandwidth (The Silent Killer)
*   **Scenario**: 13,000 VMs start simultaneously.
*   **Startup Storm**: 13,000 $\times$ Docker Pull (5GB) + Input Download. This will throttle the GCS bucket or Artifact Registry.
*   **Checkpoint Storm**: If using Option 2 (Checkpointing), 13,000 VMs writing 50MB files every 10 minutes.
    *   **Real Risk**: bandwidth is manageable (1 TB/min), but **API Rate Limiting (429 Too Many Requests)** is the killer. 13,000 synchronous `gsutil cp` calls will trigger GCS denial-of-service protection.
    *   **Mitigation**: **Jitter** (Randomized delays) is mandatory.

---

## 2. Option Analysis at Scale

### Option 1: Task Retries (No Checkpoint)
**"The Brute Force Approach"**

*   **Scaling Behavior**:
    *   **Bandwidth**: Low. Only reads at start, writes at end.
    *   **Spot Risk**: **Severe**. With 13,000 tasks, Preemption Rate is a probability game. If the rate is 5% per hour, probability of a 7-hour task surviving is $(0.95)^7 \approx 70\%$. Approximately 30% of tasks will require retries.
    *   **Cost Impact**: Costs are incurred for the 30% waste.
*   **Verdict**: Viable only if concurrency is limited (e.g., "Rolling Batch" of 50 molecules at a time) to avoid triggering aggressive preemption logic in the zone.

### Option 2: Checkpoint & Resume
**"The Robust Approach"**

*   **Scaling Behavior**:
    *   **Bandwidth**: High. The "Watchdog" syncer adds constant background noise.
    *   **Spot Risk**: **Negligible**. Preemptions are just speed bumps.
    *   **Throughput**: Higher effective throughput because progress is not lost.
*   **Implementation Detail**:
    *   **Watchdog**: A lightweight Bash background loop. No need for sidecar containers.
    *   **Jitter**: `sleep $(( 600 + RANDOM % 120 ))` is mandatory to spread API load.

---

## 3. Recommended Strategy: "Native Batch Throttling"

Do not build a custom orchestrator. Use Google Cloud Batch's built-in queue management.

### Architecture
1.  **Submission**: Submit larger jobs (e.g., 50 molecules = 650 tasks) or one massive Array Job.
2.  **Throttling**: Use the `parallelism` field in the Batch Job Spec.
    *   `taskCount`: 13,000 (Total work)
    *   `parallelism`: 500 (Max concurrent GPUs)
3.  **Result**:
    *   Google queues the 12,500 pending tasks automatically.
    *   As one task finishes, the Next task starts.
    *   **Benefit**: Keeps the system within Quota limits and smooths out the "Startup Storm" over time without writing any custom queue logic.

### The "Hybrid" Checkpoint
For the 1,000 molecule scale, **Option 2 (Checkpointing) is superior**.
*   **Reason**: At 13,000 tasks, maintenance windows, stockouts, and preemptions are inevitable. Relying on luck (Option 1) for 13,000 dice rolls is statistically poor.
*   **Cost**: The saved compute from not re-running 7-hour jobs pays for the implementation effort.
