import pandas as pd
import hashlib
import glob
import os
import sys
from datetime import datetime

# --- CONFIG ---
# Cloud Batch Mount Path
MOUNT_BASE = "/mnt/disks/gcs/output" 
OUTPUT_BASE = "/mnt/disks/gcs/output/consolidated_bq"

# Input Patterns (File System)
MANIFEST_PATH = f"{MOUNT_BASE}/consolidated/manifest.csv"
RDKIT_PATTERN = f"{MOUNT_BASE}/prod_100k/rdkit/shard_*.csv"
TXGEMMA_PATTERN = f"{MOUNT_BASE}/prod_100k/txgemma/shard_*.csv"
GNINA_PATTERN = f"{MOUNT_BASE}/prod_100k/gnina/shard_*/predictions_*.csv"

RUN_ID = "prod_100k_v1"

def load_csvs(pattern):
    """Read multiple CSVs from Filesystem into a single DataFrame."""
    files = glob.glob(pattern)
    print(f"Found {len(files)} files for pattern {pattern}")
    
    dfs = []
    for f in files:
        try:
            df = pd.read_csv(f)
            dfs.append(df)
        except Exception as e:
            print(f"Error reading {f}: {e}")
            
    if not dfs:
        return pd.DataFrame()
    
    full_df = pd.concat(dfs, ignore_index=True)
    return full_df

def generate_hash(smiles):
    """SHA256 hash of SMILES string."""
    try:
        return hashlib.sha256(str(smiles).encode()).hexdigest()
    except:
        return "INVALID_HASH"

def main():
    print("--- STARTING CONSOLIDATION (BATCH MODE) ---", flush=True)
    
    # Verify Mounts
    if not os.path.exists(MOUNT_BASE):
        print(f"CRITICAL: Mount path {MOUNT_BASE} does not exist.", flush=True)
        sys.exit(1)
        
    os.makedirs(OUTPUT_BASE, exist_ok=True)
    
    print("--- 1. Loading Manifest ---", flush=True)
    if not os.path.exists(MANIFEST_PATH):
         print(f"CRITICAL: Manifest not found at {MANIFEST_PATH}", flush=True)
         sys.exit(1)
         
    df_manifest = pd.read_csv(MANIFEST_PATH)
    print(f"Manifest: {len(df_manifest)} rows", flush=True)
    
    print("--- 2. Loading Scores ---", flush=True)
    df_rdkit = load_csvs(RDKIT_PATTERN) 
    print(f"RDKit: {len(df_rdkit)} rows", flush=True)
    
    df_txgemma = load_csvs(TXGEMMA_PATTERN) 
    print(f"TxGemma: {len(df_txgemma)} rows", flush=True)
    
    # Gnina check (Recursive glob might be needed depending on python version, 
    # but standard glob usually doesn't do ** by default in older versions. 
    # We'll stick to the specific pattern or assume shallow structure if possible.
    # Actually, Gnina uses shard_N/predictions_N.csv. Standard glob handles wildcard dirs.)
    df_gnina = load_csvs(GNINA_PATTERN) 
    print(f"Gnina: {len(df_gnina)} rows", flush=True)
    
    print("--- 3. Merging ---", flush=True)
    # Left Join everything to Manifest
    merged = df_manifest.merge(df_rdkit, on="smiles", how="left")
    merged = merged.merge(df_txgemma, on="smiles", how="left")
    merged = merged.merge(df_gnina, on="smiles", how="left")
    
    print(f"Merged Shape: {merged.shape}", flush=True)
    
    print("--- 4. Formatting for BigQuery ---", flush=True)
    merged['molecule_hash'] = merged['smiles'].apply(generate_hash)
    merged['run_id'] = RUN_ID
    merged['created_at'] = datetime.utcnow().isoformat()
    
    # --- TABLE 1: molecule_registry ---
    registry_cols = ['molecule_hash', 'smiles', 'global_id', 'sdf_gcs_path', 'run_id', 'created_at']
    # Ensure cols exist
    for c in registry_cols:
        if c not in merged.columns: merged[c] = None
        
    df_registry = merged[registry_cols].copy()
    reg_csv = f"{OUTPUT_BASE}/bq_molecule_registry.csv"
    df_registry.to_csv(reg_csv, index=False)
    print(f"Saved {reg_csv}", flush=True)
    
    # --- TABLE 2: screening_results ---
    rename_map = {
        "qed": "qed_score",
        "toxicity": "toxicity_label"
    }
    merged.rename(columns=rename_map, inplace=True)
    
    screening_cols = ['molecule_hash', 'docking_score', 'qed_score', 'sa_score', 'logp', 'mw', 'toxicity_label']
    for c in screening_cols:
        if c not in merged.columns: merged[c] = None
        
    df_screening = merged[screening_cols].copy()
    screen_csv = f"{OUTPUT_BASE}/bq_screening_results.csv"
    df_screening.to_csv(screen_csv, index=False)
    print(f"Saved {screen_csv}", flush=True)

    print("--- Done. ---", flush=True)

if __name__ == "__main__":
    main()
