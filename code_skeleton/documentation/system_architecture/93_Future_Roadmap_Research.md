# Future Roadmap: Advancing the Architecture (Research Integration)

Based on a comprehensive review of the current state of Deep Learning in Drug Discovery (referencing *papers-for-molecular-design-using-DL*), we have identified three high-impact integrations to evolve the platform from "MVP" to "State of the Art".

## 1. Upgrade the Generator: Diffusion Models
**Current:** Pocket2Mol (Autoregressive GNN).
**Upgrade:** **TargetDiff** (3D Equivariant Diffusion).

*   **The Logic:** Autoregressive models build molecules atom-by-atom, which can lead to local optima and invalid valencies. Diffusion models generate the entire molecular graph simultaneously by denoising from a prior distribution.
*   **Impact:** Higher valid generation rate, better handling of complex ring systems, and improved binding affinity distributions.
*   **Integration:** Drop-in replacement for the `Generation` component (Same Input: PDB, Same Output: SDF).

## 2. Upgrade the Critic: 3D Representation Learning
**Current:** XGBoost on 2D Morgan Fingerprints.
**Upgrade:** **Uni-Mol** (3D Molecular Transformers).

*   **The Logic:** Our current pipeline generates 3D molecules, but the "Critic" (Proxy Model) flattens them into 2D fingerprints to predict scores. This discards critical shape and electrostatic information.
*   **Impact:** By using Uni-Mol to generate **3D Embeddings** of the protein-ligand complex, the Proxy Model can "see" the fit, drastically improving the correlation between Predicted $\Delta G$ and Oracle $\Delta G$.
*   **Integration:** Replace `rdkit.GetMorganFingerprint` in the Training Component with the Uni-Mol encoder.

## 3. Upgrade the Loop: Reinforcement Learning Frameworks
**Current:** Custom "Weighted Retraining" scripts.
**Upgrade:** **REINVENT 4** Core Logic.

*   **The Logic:** Instead of maintaining custom selection/ranking scripts, adopt the reward shaping and curriculum learning algorithms standardized in REINVENT 4.
*   **Impact:** Built-in diversity filters (Scaffold Memory), simpler configuration of Multi-Parameter Optimization (MPO), and more stable policy gradient updates.
*   **Integration:** Adapt the `Selector` component to use REINVENT's scoring logic.

---

## Strategic Roadmap

| Phase | Component | Technology | Goal |
| :--- | :--- | :--- | :--- |
| **Phase 4 (Next)** | Critic | **Uni-Mol** | Fix the "2D blindspot" in the 3D generation loop. |
| **Phase 5** | Generator | **TargetDiff** | benchmark Diffusion vs Autoregressive generation. |
| **Phase 6** | Oracle | **MM-GBSA** | Replace Mock Oracle with physics-based "Light Oracle" (GROMACS). |
