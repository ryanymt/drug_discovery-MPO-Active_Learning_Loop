import os
import csv
import subprocess
import re
import argparse
import sys

# Default Configuration
DEFAULT_GCS_BASE = "gs://ryanymt/output/generated"
DEFAULT_OUTPUT_FILE = "manifest.csv"
DEFAULT_GCS_OUTPUT_PATH = "gs://ryanymt/output/consolidated/manifest.csv"

def get_shard_files(gcs_base):
    """List all SMILES.txt files using gsutil"""
    print(f"Scanning {gcs_base} for SMILES.txt files...")
    # Use gsutil to list all SMILES.txt files
    cmd = f"gsutil ls '{gcs_base}/shard_*/SMILES.txt'"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"Error listing files: {result.stderr}")
        return []
    
    files = result.stdout.strip().split('\n')
    # Filter empty strings if any
    files = [f for f in files if f.strip()]
    print(f"Found {len(files)} shards.")
    return files

def process_shard(gcs_path, gcs_base):
    """Read a single SMILES.txt and return list of records"""
    # Extract shard ID from path: .../shard_12/SMILES.txt
    match = re.search(r'shard_(\d+)', gcs_path)
    if not match:
        print(f"Could not parse shard ID from {gcs_path}")
        return []
    
    shard_id = match.group(1)
    
    # Read content
    cmd = f"gsutil cat '{gcs_path}'"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"Error reading {gcs_path}: {result.stderr}")
        return []
    
    lines = result.stdout.strip().split('\n')
    records = []
    
    for i, smile in enumerate(lines):
        if not smile.strip():
            continue
            
        local_index = i
        # Use a cycle-aware prefix if possible, or just standard gen_shard
        global_id = f"gen_shard{shard_id}_idx{local_index}"
        sdf_path = f"{gcs_base}/shard_{shard_id}/SDF/{local_index}.sdf"
        
        records.append({
            "global_id": global_id,
            "smiles": smile.strip(),
            "shard_id": shard_id,
            "local_index": local_index,
            "sdf_gcs_path": sdf_path
        })
        
    return records

def main():
    parser = argparse.ArgumentParser(description="Consolidate distributed shard manifests.")
    parser.add_argument("--input_base", default=DEFAULT_GCS_BASE, help="GCS Base path (e.g., gs://bucket/output/generated)")
    parser.add_argument("--output_csv", default=DEFAULT_OUTPUT_FILE, help="Local output CSV filename")
    parser.add_argument("--upload_to", default=DEFAULT_GCS_OUTPUT_PATH, help="GCS Destination for consolidated CSV")
    
    args = parser.parse_args()
    
    files = get_shard_files(args.input_base)
    if not files:
        print("No files found. Exiting.")
        return

    all_records = []
    
    print(f"Processing {len(files)} shards sequentially...")
    
    for i, file in enumerate(files):
        if i % 10 == 0:
            print(f"Processed {i}/{len(files)} shards...")
            
        records = process_shard(file, args.input_base)
        all_records.extend(records)

    # Sort logic
    all_records.sort(key=lambda x: (int(x['shard_id']), int(x['local_index'])))
    
    print(f"Total molecules consolidated: {len(all_records)}")
    
    # Write to CSV
    print(f"Writing to {args.output_csv}...")
    with open(args.output_csv, 'w', newline='') as f:
        fieldnames = ['global_id', 'smiles', 'shard_id', 'local_index', 'sdf_gcs_path']
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(all_records)
        
    # Upload to GCS
    print(f"Uploading to {args.upload_to}...")
    subprocess.run(f"gsutil cp {args.output_csv} {args.upload_to}", shell=True, check=True)
    print("Done.")

if __name__ == "__main__":
    main()
