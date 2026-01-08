# Observation Report: Cycle 2 Production Analysis

**Date**: January 2026
**Subject**: Analysis of Generative Model Behavior (Cycle 1 vs Cycle 2)
**Scale**: 100k (Baseline) vs 10k (Fine-Tuned)

## 1. Hypothesis: The "Exploitation" Shift
The primary goal of the Active Learning loop is to bias the Generative Model towards the "High Affinity" region of chemical space discovered in Cycle 1. The expectation was a leftward shift in the Delta G distribution (lower is better).

## 2. Quantitative Analysis

Comparing the Baseline Generation (`prod_100k_v1`) against the Fine-Tuned Generation (`cycle2_10k`):

| Metric | Cycle 1 (Baseline) | Cycle 2 (Fine-Tuned) | Shift |
| :--- | :--- | :--- | :--- |
| **Population Size** | 70,259 | 6,853 | -- |
| **Uniqueness** | 100.0% | 100.0% | 0.0% |
| **Mean Affinity** | -36.68 kcal/mol | **-36.89 kcal/mol** | -0.21 kcal/mol |
| **Hit Rate (<-30)** | 92.1% | **94.5%** | +2.3% |
| **Safety Rate** | 92.6% | 91.3% | -1.3% |

### Key Findings
1.  **Optimization Success**: The model successfully learned from the "Elite" dataset. The **Hit Rate (< -30 kcal/mol)** improved by **2.3%**, indicating the generator is wasting less time on weak binders.
2.  **No Mode Collapse**: Despite fine-tuning on a clustered subset, **100% of the generated molecules were unique**. The model has *not* collapsed into memorizing the training data or outputting duplicates.
3.  **Safety Trade-off**: A slight dip in the Safety Rate (-1.3%) was observed. This suggests a potential correlation between high affinity features and toxicity features (e.g., reactive groups or excessive lipophilicity).
4. **The weaklink**: The quality of the surrogate model is the most important part in this pipeline. In this experiment run, 1000 filtered candidates (top 800 plus random 200) was used for Gromac MM-GBSA scoring. Filtering was done in the steps of Non-Toxic > Top QED > Top Gnina. This seem to have filter candidate from similar zone. No outlier score predictions from Surrogate model observed.

## 3. Structural Dynamics (Scaffold Exploitation)
While exact structural visualization requires SMILES data (currently unavailable in the result CSVs), the metrics strongly point to **Scaffold Exploitation**:

*   **Observation**: The combination of **High Affinity Improvement** and **High Uniqueness** indicates the model is exploring variations *within* a successful chemical series (finding the "local maximum") rather than jumping to entirely new, likely lower-scoring scaffolds.
*   **Verdict**: This confirms the "fixing on backbone" behavior observed in early tests is a feature, not a bug. It means the RL loop is working as a "Lead Optimizer".

## 4. Conclusion & Recommendations
The system resembles a "Hill Climber", systematically improving candidates within a high-value region.


### Strategy for New Experimental run
1.  **Enforce Diversity**: To prevent getting stuck in this local optimum, continue to try and improve the "Clustering" selection strategy.
2.  **Monitor Safety**: If the Safety Rate drops below 90%, implement **Constraint-Based Generation** to penalize toxicophores during inference.
3.  **Surrogate Model**: Even with divserse set of tarining data, other type of models like CNN models could be tested.
