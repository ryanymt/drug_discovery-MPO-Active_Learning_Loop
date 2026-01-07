import subprocess
import csv
import os
import glob
import math

# --- Configuration ---
PROJECT_ID = "gcda-apac-sc"
BQ_SOURCE_TABLE = "bioops_platform.elite_candidates_v1" # The 10k Elite
OUTPUT_BUCKET_BASE = "gs://ryanymt/output/generated"
DEST_BUCKET = "gs://ryanymt/input/training/elite_sdfs"
LOCAL_WORK_DIR = "elite_prep"

def run_command(cmd, shell=False):
    """Runs a shell command and returns output."""
    try:
        if shell:
            result = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
        else:
            result = subprocess.check_output(cmd, stderr=subprocess.STDOUT)
        return result.decode("utf-8").strip()
    except subprocess.CalledProcessError as e:
        # print(f"Error running command: {cmd}") 
        # Output might be too verbose to print always
        return None

def get_candidates():
    """Fetches candidates from BigQuery (Elite JOIN Registry)."""
    print("Fetching elite candidates from BigQuery...")
    # We need global_id to find the physical file (shard_id / index)
    query = f"""
        SELECT DISTINCT t1.molecule_hash, t2.global_id 
        FROM `{PROJECT_ID}.{BQ_SOURCE_TABLE}` t1
        JOIN `{PROJECT_ID}.bioops_platform.molecule_registry` t2
        ON t1.molecule_hash = t2.molecule_hash
    """
    cmd = ["bq", "query", "--quiet", "--format=csv", "--nouse_legacy_sql", "--max_rows=15000", query]
    output = run_command(cmd)
    
    candidates = []
    if output:
        lines = output.splitlines()
        reader = csv.DictReader(lines)
        for row in reader:
            # Parse global_id: gen_shard0_idx1
            gid = row['global_id']
            try:
                parts = gid.split('_')
                shard = parts[1].replace('shard', '')
                idx = parts[2].replace('idx', '')
                row['shard_id'] = shard
                row['local_index'] = idx
                candidates.append(row)
            except:
                print(f"Failed to parse global_id: {gid}")
    return candidates

def find_sdf_base_path(shard_id):
    """Finds the actual SDF directory for a given shard."""
    base_search = f"{OUTPUT_BUCKET_BASE}/shard_{shard_id}/"
    ls_cmd = f"gsutil ls {base_search}"
    output = run_command(ls_cmd, shell=True)
    
    if not output:
        return None
        
    paths = output.splitlines()
    # Look for run_config directories
    run_configs = [p for p in paths if "run_config" in p]
    run_configs.sort()
    
    for rc in reversed(run_configs):
        sdf_path = f"{rc}SDF/"
        check_cmd = f"gsutil ls {sdf_path}"
        if run_command(check_cmd, shell=True):
            return sdf_path
            
    return None

def main():
    os.makedirs(LOCAL_WORK_DIR, exist_ok=True)
    
    # 1. Get Candidates
    candidates = get_candidates()
    print(f"Found {len(candidates)} candidates.")
    
    if len(candidates) == 0:
        print("No candidates found via BQ. Exiting.")
        return

    # 2. Map Shards to SDF Paths
    shard_map = {}
    unique_shards = set(c['shard_id'] for c in candidates)
    print(f"Scanning {len(unique_shards)} unique shards for SDF paths...")
    
    for s_id in unique_shards:
        path = find_sdf_base_path(s_id)
        if path:
            shard_map[s_id] = path
            # print(f"Shard {s_id} -> {path}")
        else:
            print(f"WARNING: Could not find SDF path for Shard {s_id}")

    # 3. Download and Aggregate
    # We will simply download ALL valid SDFs into flat structure or batches.
    # Pocket2Mol Training usually processes raw SDFs. 
    # Let's organize by batch to avoid directory overload? 
    # Actually, for training dataset creation, we might want to consolidate them.
    # Let's save them as `hash.sdf` in the destination bucket.
    # BUT 10k files is slow for GSUtil. 
    # Better: Keep batch structure (e.g. 100 SDfs per file) and let the training script unroll them?
    # Or just download locally here and create a TAR?
    # Let's try parallel download to local dir, then TAR.
    
    print("Downloading SDFs...")
    success_count = 0
    
    # Group by shard to minimize 'gsutil cp' calls (bulk copy?)
    # Can't bulk copy easily because filenames are 'idx.sdf' and we want 'hash.sdf' or similar?
    # Actually, 'idx.sdf' is fine if we keep them separate.
    # Let's iterate and download.
    
    # Optimization: Generate a list of (src, dst) and run parallel gsutil?
    # Or simpler: Just download sequentially/batches for now. 10k is manageable.
    
    # Optimized Parallel Download using Python Client
    print("Initializing GCS Client...")
    try:
        from google.cloud import storage
        storage_client = storage.Client()
    except ImportError:
        print("Error: google-cloud-storage not installed. Please install it.")
        return

    tasks = []
    print("Collecting download tasks...")
    
    for cand in candidates:
        shard = cand['shard_id']
        idx = cand['local_index']
        mol_hash = cand['molecule_hash']
        
        if shard not in shard_map: continue
        
        # gs://bucket/path/to/file -> bucket, blob_name
        src_url = f"{shard_map[shard]}{idx}.sdf"
        
        # Parse bucket and blob
        # Expected format: gs://bucket/path/obj
        parts = src_url.replace("gs://", "").split("/", 1)
        bucket_name = parts[0]
        blob_name = parts[1]
        
        dst_local = os.path.join(LOCAL_WORK_DIR, f"{mol_hash}.sdf")
        tasks.append((bucket_name, blob_name, dst_local))

    print(f"Queueing {len(tasks)} downloads...")
    
    # Threaded Download
    from concurrent.futures import ThreadPoolExecutor
    
    success_count = 0
    
    def download_blob(task):
        b_name, b_path, dst = task
        if os.path.exists(dst): return True
        try:
            bucket = storage_client.bucket(b_name)
            blob = bucket.blob(b_path)
            blob.download_to_filename(dst)
            return True
        except Exception as e:
            print(f"Failed {b_path}: {e}")
            return False

    with ThreadPoolExecutor(max_workers=20) as executor:
        results = list(executor.map(download_blob, tasks))
    
    success_count = sum(results)
    print(f"Downloaded {success_count} SDFs.")
    
    # Cleanup main
    # ...
    
    # 4. Tar and Upload
    tar_name = "elite_sdfs.tar.gz"
    print(f"Creating archive {tar_name}...")
    run_command(f"tar -czf {tar_name} -C {LOCAL_WORK_DIR} .", shell=True)
    
    dest_url = f"{DEST_BUCKET}/{tar_name}"
    print(f"Uploading to {dest_url}...")
    run_command(f"gsutil cp {tar_name} {dest_url}", shell=True)
    
    print("Done.")

if __name__ == "__main__":
    main()
