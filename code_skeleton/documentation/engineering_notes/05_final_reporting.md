# 05. Final Reporting & Visualization

The "Hero Plots" demonstrating the success of the Active Learning loops are generated using a standalone Python script located in the workspace.

## Script: `workspace/plot_pareto.py`

This script compares two datasets (e.g., "Baseline/Cycle 1" vs "Active Learning/Cycle N") and generates statistical visualizations showing the improvement in Binding Affinity and Drug-likeness.

### Prerequisites
*   **Environment**: Python 3.x
*   **Dependencies**: `pandas`, `matplotlib`, `seaborn`.
    *   *Note*: The script attempts to auto-install these if missing.
*   **Inputs**: Two CSV files containing minimal columns:
    *   `gnina_cnn_affinity`: Binding prediction (pK/Affinity).
    *   `rdkit_qed`: Drug-likeness score.
    *   `rdkit_sa_score`: Synthesizability score.

### Usage
Run the script locally or on a cloud VM after downloading the results from GCS.

```bash
# 1. Download Results (Example)
gsutil cp gs://drug-discovery-mvp-docking-results/output/loop_1/selection/joined_results.csv baseline.csv
gsutil cp gs://drug-discovery-mvp-docking-results/output/loop_2/selection/joined_results.csv loop2.csv

# 2. Run Plotter
python3 workspace/plot_pareto.py \
    --baseline baseline.csv \
    --loop2 loop2.csv \
    --output ./report_plots
```

### Outputs
The script generates three key images in the output directory:

| Filename | Description | Purpose |
| :--- | :--- | :--- |
| `scatter_affinity_qed.png` | **Optimization Landscape** <br> (Scatter Plot) | Shows the "Pareto Frontier" expansion. Look for new points in the top-right quadrant (High Affinity + High QED). |
| `kde_affinity_shift.png` | **Distribution Shift** <br> (Density Plot) | Demonstrates the "Hill Climbing" effect. The Loop 2 curve should be shifted to the right (higher affinity) compared to Baseline. |
| `box_metrics_comparison.png` | **Statistical Summary** <br> (Box Plots) | Side-by-side comparison of median scores for Affinity, QED, and SA. |
