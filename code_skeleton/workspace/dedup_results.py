
import csv

input_file = "bq_oracle_results.csv"
output_file = "bq_oracle_results_clean.csv"

# Dict to store best score per hash
# key: molecule_hash, value: (delta_g, row_dict)
best_results = {}

print("Deduping results...")

with open(input_file, 'r') as f:
    reader = csv.DictReader(f)
    fieldnames = reader.fieldnames
    
    for row in reader:
        h = row['molecule_hash']
        delta_g_str = row['delta_g']
        
        # Skip failed rows if we have successes? Or keep them?
        # Better to keep successes.
        if row['status'] != 'Success' or delta_g_str == 'NaN':
            # Store failure only if we don't have a success yet
            if h not in best_results:
                best_results[h] = (9999.0, row) 
            continue
            
        try:
            delta_g = float(delta_g_str)
        except:
            continue
            
        if h not in best_results:
            best_results[h] = (delta_g, row)
        else:
            # Update if better (lower) score
            if delta_g < best_results[h][0]:
                best_results[h] = (delta_g, row)

# Write output
with open(output_file, 'w') as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    for _, (_, row) in best_results.items():
        writer.writerow(row)

print(f"Original rows: Unknown, Unique Hashes: {len(best_results)}")
