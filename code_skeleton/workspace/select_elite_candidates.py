import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
from rdkit.ML.Cluster import Butina
import time
import argparse
import sys
import os

# Parameters
TARGET_CLUSTERS = 8000
TARGET_RANDOM = 2000
# Tighten cutoff to 0.35 (sim > 0.65) to get more clusters (smaller clusters)
CLUSTERING_CUTOFF = 0.35 

def cluster_fingerprints(fps, cutoff=0.5):
    nfps = len(fps)
    print(f"Calculating distances for {nfps} fingerprints...")
    dists = []
    # Bulk Tanimoto is optimized in C++
    for i in range(1, nfps):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])
        dists.extend([1-x for x in sims])
    
    print("Clustering (Butina)...")
    start = time.time()
    clusters = Butina.ClusterData(dists, nfps, distThresh=cutoff, isDistData=True, reordering=True)
    end = time.time()
    print(f"Clustering complete in {end-start:.2f}s. Found {len(clusters)} clusters.")
    return clusters

def select_elite(df):
    # 1. Ranking Pool
    # Already sorted by final_score ASC in input
    pool = df.copy()
    pool.reset_index(drop=True, inplace=True)
    
    print(f"Generating fingerprints for pool of {len(pool)}...")
    mols = [Chem.MolFromSmiles(s) for s in pool.smiles]
    fps = [AllChem.GetMorganFingerprintAsBitVect(x, 2, 1024) for x in mols]
    
    # 2. Clustering
    clusters = cluster_fingerprints(fps, cutoff=CLUSTERING_CUTOFF)
    
    # 3. Select Best from Clusters
    selection_indices = []
    print(f"Selecting best candidates from clusters (Target: {TARGET_CLUSTERS})...")
    
    # For each cluster (list of indices from pool), pick the one with MIN index.
    # Since pool is sorted by score ASC, min index = best score.
    cluster_reps = []
    for cluster in clusters:
        best_idx = min(cluster)
        cluster_reps.append(best_idx)
    
    cluster_reps.sort()
    
    if len(cluster_reps) > TARGET_CLUSTERS:
        selected_indices = cluster_reps[:TARGET_CLUSTERS]
        print(f"Selected {len(selected_indices)} structured candidates (Hit Cap).")
    else:
        selected_indices = cluster_reps
        print(f"Selected {len(selected_indices)} structured candidates (All Clusters).")
        
    actual_structured_count = len(selected_indices)
    shortfall = TARGET_CLUSTERS - actual_structured_count
    
    if shortfall > 0:
        print(f"Shortfall of {shortfall} candidates from clustering. Will backfill with diversity.")
        # We will increase the diversity target dynamically
        global TARGET_RANDOM
        TARGET_RANDOM += shortfall
    
    print(f"Selected {len(selected_indices)} structured candidates.")
    
    elite_pool = pool.iloc[selected_indices].copy()
    elite_pool['selection_strategy'] = 'cluster_centroid'
    
    # 4. Diversity Mix (Random 2000 from remainder)
    # Since we only loaded 20k into 'pool', we don't have the full remainder here.
    # But wait, User said: "Add 2,000 random molecules."
    # My BQ query only exported 20,000.
    # So I have to select randoms from the *remaining* 12,000 of this pool?
    # OR did user mean random from the WHOLE filtered set?
    # User: "The Diversity Mix: Add your 2,000 random molecules." context: "clustering these 20,000 ... pick best ... add 2000".
    # It implies adding to the selection. If the 20k was the "Survivors check", then sure.
    # NOTE: To keep it simple and safe within this job, we will pick randoms from the UNSELECTED portion of the 20k pool.
    # This 20k pool represents the top tier. Random sampling within top tier is "Local Diversity".
    
    selected_indices_set = set(selected_indices)
    remaining_indices = [i for i in range(len(pool)) if i not in selected_indices_set]
    
    # If we need more diversity, we might need to query BQ again, but script is standalone.
    # Let's sample from the remainder of the 20k.
    
    remainder_df = pool.iloc[remaining_indices]
    print(f"Selecting {TARGET_RANDOM} random diversity candidates from remainder of top 20k ({len(remainder_df)})...")
    
    if len(remainder_df) > TARGET_RANDOM:
        random_picks = remainder_df.sample(TARGET_RANDOM)
    else:
        random_picks = remainder_df
        
    random_picks = random_picks.copy()
    random_picks['selection_strategy'] = 'diversity_random'
    
    # 5. Combine
    final_elite = pd.concat([elite_pool, random_picks])
    return final_elite

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    
    print(f"Reading input from {args.input}...")
    df = pd.read_csv(args.input)
    
    # Select
    elite_df = select_elite(df)
    
    # Export
    print(f"Exporting {len(elite_df)} rows to {args.output}...")
    
    # Ensure parent dir exists (Safety for FUSE/Local)
    output_dir = os.path.dirname(args.output)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
        
    elite_df.to_csv(args.output, index=False)
    print("Done.")

if __name__ == "__main__":
    main()
