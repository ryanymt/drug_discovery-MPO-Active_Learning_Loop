import subprocess
import csv
import os
import glob
import math

# --- Configuration ---
PROJECT_ID = "gcda-apac-sc"
BQ_TABLE = "bioops_platform.oracle_candidates_v1"
OUTPUT_BUCKET_BASE = "gs://ryanymt/output/generated"
DEST_BUCKET = "gs://ryanymt/input/oracle_batches"
LOCAL_WORK_DIR = "oracle_prep"

def run_command(cmd, shell=False):
    """Runs a shell command and returns output."""
    try:
        if shell:
            result = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
        else:
            result = subprocess.check_output(cmd, stderr=subprocess.STDOUT)
        return result.decode("utf-8").strip()
    except subprocess.CalledProcessError as e:
        print(f"Error running command: {cmd}")
        print(e.output.decode("utf-8"))
        return None

def get_candidates():
    """Fetches candidates from BigQuery."""
    print("Fetching candidates from BigQuery...")
    query = f"""
        SELECT DISTINCT t1.molecule_hash, t2.global_id 
        FROM `{PROJECT_ID}.{BQ_TABLE}` t1
        JOIN `{PROJECT_ID}.bioops_platform.molecule_registry` t2
        ON t1.molecule_hash = t2.molecule_hash
    """
    cmd = ["bq", "query", "--quiet", "--format=csv", "--nouse_legacy_sql", "--max_rows=2000", query]
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
    # List all directories in the shard
    ls_cmd = f"gsutil ls {base_search}"
    output = run_command(ls_cmd, shell=True)
    
    if not output:
        return None
        
    paths = output.splitlines()
    # Look for run_config directories
    run_configs = [p for p in paths if "run_config" in p]
    
    # Sort to hopefully get the latest? Or check which one has SDFs.
    # Usually the latest timestamp is last.
    run_configs.sort()
    
    for rc in reversed(run_configs):
        # Check if SDF dir exists
        sdf_path = f"{rc}SDF/"
        # We can try to list one file to verify, or just trust the structure if directory exists
        # Let's try to ls the directory
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
        print("No candidates found. Exiting.")
        return

    # 2. Map Shards to SDF Paths
    shard_map = {}
    unique_shards = set(c['shard_id'] for c in candidates)
    print(f"Scanning {len(unique_shards)} unique shards for SDF paths...")
    
    for s_id in unique_shards:
        path = find_sdf_base_path(s_id)
        if path:
            shard_map[s_id] = path
            print(f"Shard {s_id} -> {path}")
        else:
            print(f"WARNING: Could not find SDF path for Shard {s_id}")

    # 3. Download and Aggregate
    # User Request: 1000 tasks, 1 mol/task.
    batch_size = 1
    total_batches = math.ceil(len(candidates) / batch_size)
    
    print(f"Preparing {total_batches} batches...")
    
    mapping_rows = []

    for i in range(total_batches):
        batch_candidates = candidates[i*batch_size : (i+1)*batch_size]
        batch_filename = f"batch_{i}.sdf"
        local_batch_path = os.path.join(LOCAL_WORK_DIR, batch_filename)
        
        print(f"Processing Batch {i}: {len(batch_candidates)} molecules...")
        
        with open(local_batch_path, 'w') as batch_file:
            for b_idx, cand in enumerate(batch_candidates):
                shard = cand['shard_id']
                idx = cand['local_index']
                mol_hash = cand['molecule_hash']
                
                # Add to mapping
                # Batch ID (i), Mol Index in Batch (b_idx), Hash
                mapping_rows.append({"batch_id": i, "mol_index": f"mol_{b_idx}", "molecule_hash": mol_hash})
                
                if shard not in shard_map:
                    print(f"Skipping mol {shard}/{idx} - Path not found.")
                    continue
                    
                sdf_url = f"{shard_map[shard]}{idx}.sdf"
                local_temp = os.path.join(LOCAL_WORK_DIR, f"temp_{shard}_{idx}.sdf")
                
                # Download
                cp_cmd = f"gsutil cp {sdf_url} {local_temp}"
                if run_command(cp_cmd, shell=True):
                    # Append to batch file
                    with open(local_temp, 'r') as f:
                        batch_file.write(f.read())
                        batch_file.write("\n$$$$\n") 
                    os.remove(local_temp)
                else:
                    print(f"Failed to download {sdf_url}")
                    
        # Upload Batch
        dest_url = f"{DEST_BUCKET}/{batch_filename}"
        run_command(f"gsutil cp {local_batch_path} {dest_url}", shell=True)
        print(f"Uploaded {dest_url}")
        os.remove(local_batch_path)

    # 4. Write and Upload Manifest
    manifest_path = os.path.join(LOCAL_WORK_DIR, "oracle_batch_map.csv")
    with open(manifest_path, 'w') as f:
        writer = csv.DictWriter(f, fieldnames=["batch_id", "mol_index", "molecule_hash"])
        writer.writeheader()
        writer.writerows(mapping_rows)
    
    manifest_dest = f"{DEST_BUCKET}/oracle_batch_map.csv"
    run_command(f"gsutil cp {manifest_path} {manifest_dest}", shell=True)
    print(f"Uploaded Mapping: {manifest_dest}")
    os.remove(manifest_path)

if __name__ == "__main__":
    main()
