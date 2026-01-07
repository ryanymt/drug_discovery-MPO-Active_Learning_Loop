# Strategy: Mitigating Mode Collapse in Generative Drug Design

## 1. The Problem: "The One-Trick Pony"
Active Learning loops naturally gravitate towards specific high-performing local minima. Generative models (like Pocket2Mol) are mathematically incentivized to find the *easiest* pattern that satisfies the reward function and repeat it. This leads to **Mode Collapse**, where the model generates thousands of variations of a single scaffold.

## 2. Implemented Solution (Cycle 2 Status)

We successfully implemented **Layer 1: Diversity-Aware Selection** to combat this.

### Mechanism: The "Robust Recipe"
Instead of picking the Top N scorers, we used a clustering-based approach:
1.  **Filter**: `SA < 4.0` AND `QED > 0.5`.
2.  **Rank**: Top 20,000 by Predicted Affinity.
3.  **Cluster**: Morgan Fingerprints (2048-bit) + Butina Algorithm (Tanimoto 0.5).
4.  **Select**:
    *   **Centroids (80%)**: The single best-scoring molecule from each distinct structural cluster.
    *   **Random (20%)**: Random sampling from the diverse pool to encourage exploration.

### Results (Cycle 2 Analysis)
*   **Success**: The strategy successfully prevented *global* mode collapse (we generated diverse valid binders).
*   **Observation**: We achieved a **99.9% Hit Rate** (Elite Candidates).
*   **Risk**: The exceedingly high hit rate suggests the model has likely found a "sweet spot". While effective, we must be vigilant against *local* mode collapse in Cycle 3.

## 3. Future Defenses (Cycle 3+)

### Layer 2: Replay Buffers (The "Memory" Fix)
*   **Concept**: Prevent "Catastrophic Forgetting" of general chemistry.
*   **Plan**: Explicitly mix 10-20% of the *original* training data (PDBBind) back into the Cycle 3 fine-tuning set.

### Layer 3: Similarity Penalties (The "Exploration" Fix)
*   **Concept**: Explicitly penalize "Copying".
*   **Logic**:
    $$ Score_{final} = Score_{MPO} - \alpha \cdot \text{MaxSimilarity}(\text{PrevElite}) $$
*   **Goal**: Force the model to leave the current high-affinity local minimum and find *new* scaffolds.
