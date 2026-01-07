
import csv
import glob
import os

def load_batch_map(map_file):
    """
    Loads the oracle_batch_map.csv into a dictionary.
    Format: batch_id,mol_index,molecule_hash
    Returns: dict { 'batch_id': 'molecule_hash' }
    """
    batch_map = {}
    with open(map_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            # We assume 1 molecule per batch for this run strategy
            b_id = row['batch_id']
            # Remove 'batch_' prefix if present in the map but not in filename logic, 
            # or handle standardized formatting.
            # The CSV likely has "batch_0", "batch_1" etc.
            batch_map[b_id] = row['molecule_hash']
    return batch_map

def consolidate_results(results_dir, batch_map, output_file):
    """
    Reads all result CSVs, maps them to hashes, and writes the final BQ import CSV.
    """
    # Pattern to match the downloaded result files
    # We downloaded them to 'results_dump/' in the previous step
    # They are named like 'batch_0_results.csv'
    file_pattern = os.path.join(results_dir, "*_results.csv")
    files = glob.glob(file_pattern)
    
    print(f"Found {len(files)} result files.")
    
    merged_data = []
    
    for filepath in files:
        filename = os.path.basename(filepath)
        # Extract batch_id from filename: batch_123_results.csv -> batch_123
        # Assuming format 'batch_N_results.csv'
        batch_id = filename.replace("_results.csv", "").replace("batch_", "")
        
        # Checking if batch_id exists in map
        if batch_id not in batch_map:
            print(f"Warning: {batch_id} not found in batch map. Skipping.")
            continue
            
        molecule_hash = batch_map[batch_id]
        
        try:
            with open(filepath, 'r') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    # Each file should have 1 row for 1-to-1 task
                    # Schema: SMILES,DeltaG,Status
                    delta_g = row.get('DeltaG', 'NaN')
                    status = row.get('Status', 'Unknown')
                    
                    # Store tuple
                    merged_data.append({
                        'molecule_hash': molecule_hash,
                        'delta_g': delta_g,
                        'status': status,
                        'batch_id': batch_id
                    })
        except Exception as e:
            print(f"Error reading {filename}: {e}")

    # Write Output
    print(f"Writing {len(merged_data)} rows to {output_file}...")
    
    headers = ['molecule_hash', 'delta_g', 'status', 'batch_id']
    with open(output_file, 'w') as f:
        writer = csv.DictWriter(f, fieldnames=headers)
        writer.writeheader()
        writer.writerows(merged_data)
        
    print("Consolidation Complete.")

if __name__ == "__main__":
    # Adjust paths as needed
    MAP_FILE = "oracle_batch_map.csv" 
    RESULTS_DIR = "../results_dump" # Relative to workspace/
    OUTPUT_FILE = "bq_oracle_results.csv"
    
    if not os.path.exists(MAP_FILE):
        print(f"Error: {MAP_FILE} not found.")
    else:
        mapping = load_batch_map(MAP_FILE)
        consolidate_results(RESULTS_DIR, mapping, OUTPUT_FILE)
