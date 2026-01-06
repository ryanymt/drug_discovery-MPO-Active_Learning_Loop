# Strategy: Mitigating Mode Collapse in Generative Drug Design

## 1. The Problem: "The One-Trick Pony"
Active Learning loops typically optimize for a score (Affinity/QED). Generative models like Pocket2Mol are mathematically incentivized to find the *easiest* pattern that satisfies this score and repeat it endlessly. This leads to **Mode Collapse**, where the model generates thousands of variations of a single scaffold, ignoring the vast diversity of chemical space.

**Why MPO isn't enough:**
Multi-Parameter Optimization (MPO) defines *what* a good molecule looks like (High Affinity + Safe), but it doesn't enforce *variety*. If one specific benzene-ring derivative satisfies all MPO constraints, the model will collapse onto it.

## 2. Solution Architecture

A three-layered defense strategy is implemented:

### Layer 1: Diversity-Aware Selection (The "Filter" Fix)
*   **Concept:** Don't just pick the Top 100 scoring molecules. Pick the **Top 100 Distinct Scaffolds**.
*   **Mechanism:**
    1.  **Cluster:** Use RDKit (Butina or MaxMin algorithm) to cluster the top 1,000 candidates by Tanimoto Similarity (Threshold 0.4).
    2.  **Representative Selection:** Pick the highest-scoring molecule from *each* cluster until the budget (100) is filled.
*   **Result:** The "Elite Set" used for training contains 100 structurally different ways to solve the binding problem.

### Layer 2: Replay Buffers (The "Memory" Fix)
*   **Concept:** Prevent "Catastrophic Forgetting" by keeping the model grounded in general chemistry.
*   **Mechanism:**
    *   **Elite Set:** The 100 new diverse winners.
    *   **Replay Set:** 20-50 molecules randomly sampled from the *original* pre-training dataset (PDBBind) or Cycle 0.
    *   **Training Mix:** Fine-tune on `Elite (80%) + Replay (20%)`.
*   **Result:** The model learns the new binding patterns without forgetting basic valency rules or chemical diversity.

### Layer 3: Similarity Penalties (The "Exploration" Fix)
*   **Concept:** Explicitly penalize "Copying".
*   **Mechanism:** During the **Selection Phase**, calculate the Tanimoto Similarity of a candidate against the *previous* cycle's Elite Set.
*   **Logic:**
    $$ Score_{final} = Score_{MPO} - \alpha \cdot \text{MaxSimilarity}(\text{PrevElite}) $$
*   **Result:** Forces the model to innovate. If it generates the same molecule as last week, it gets a low score, even if it binds well.

## 3. Implementation Roadmap

### Phase 1: Diversity Filter (Immediate)
Update the `Selector` component to implement RDKit Clustering. This is the highest ROI fix.

### Phase 2: Replay Buffer (Next Cycle)
Update the `Data Prep` script to merge a static `priors.sdf` file with the `redocked_elite.sdf` file before creating the training index.

### Phase 3: Uncertainty Sampling (Advanced)
Train an **Ensemble of 5 XGBoost Models** instead of one. Use the *variance* between their predictions as an "Uncertainty Score". Prioritize selecting molecules where the model is *confused* (High Variance) rather than just where it is confident (High Score).
