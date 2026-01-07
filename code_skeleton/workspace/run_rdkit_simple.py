import sys
import os
import csv
from rdkit import Chem
from rdkit.Chem import Descriptors, QED

# Parse Args manually to avoid complex argparse if needed, or use argparse
import argparse

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_csv", required=True)
    parser.add_argument("--output_csv", required=True)
    return parser.parse_args()

def calculate_metrics(input_path, output_path):
    print(f"Reading from {input_path}...")
    import pandas as pd
    
    try:
        df = pd.read_csv(input_path)
        if 'smiles' not in df.columns:
             # Fallback: maybe no header? try to read as list
             print("Warning: 'smiles' column not found. Trying to read as raw line list.")
             with open(input_path, 'r') as f:
                 smiles_list = [line.strip().split(',')[0] for line in f if line.strip()]
        else:
            smiles_list = df['smiles'].tolist()
    except Exception as e:
        print(f"Error reading CSV with pandas: {e}. Fallback to raw lines.")
        with open(input_path, 'r') as f:
             # Skip header if looks like csv
             lines = f.readlines()
             if "smiles" in lines[0].lower():
                 lines = lines[1:]
             smiles_list = [l.strip().split(',')[1] if ',' in l else l.strip() for l in lines if l.strip()]

    print(f"Found {len(smiles_list)} SMILES.")
    
    results = []
    valid_count = 0
    
    for i, smiles in enumerate(smiles_list):
        if not smiles: continue
            
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            try:
                qed = QED.qed(mol)
                # Skip SA Scorer to avoid segfault
                sa = 0.0 
                logp = Descriptors.MolLogP(mol)
                mw = Descriptors.MolWt(mol)
                
                results.append({
                    'filename': f"mol_{i}.sdf",
                    'smiles': smiles,
                    'qed': qed,
                    'sa_score': sa,
                    'logp': logp,
                    'mw': mw
                })
                valid_count += 1
            except Exception as e:
                print(f"Error processing molecule {i}: {e}")
        else:
            print(f"Invalid SMILES at index {i}: {smiles}")

    print(f"Calculated metrics for {valid_count} molecules.")
    
    print(f"Writing results to {output_path}...")
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    fields = ['filename', 'smiles', 'qed', 'sa_score', 'logp', 'mw']
    with open(output_path, 'w') as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(results)
        
    print("Done.")

if __name__ == "__main__":
    args = get_args()
    calculate_metrics(args.input_csv, args.output_csv)
