import argparse
import pandas as pd
import numpy as np
import xgboost as xgb
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, r2_score
import os
import joblib

# RDKit Imports
from rdkit import Chem
from rdkit.Chem import AllChem

def generate_fingerprints(smiles_list, n_bits=2048):
    """Generates ECFP4 fingerprints for a list of SMILES."""
    fps = []
    valid_indices = []
    
    print(f"Generating fingerprints for {len(smiles_list)} molecules...")
    for i, smile in enumerate(smiles_list):
        try:
            mol = Chem.MolFromSmiles(smile)
            if mol:
                fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=n_bits)
                fps.append(list(fp))
                valid_indices.append(i)
        except Exception as e:
            continue
            
    return np.array(fps), valid_indices

def train_xgboost(train_csv, output_model):
    print("--- XGBoost Proxy Model Training (Fingerprint-based) ---")
    
    # 1. Load Training Data
    print(f"Loading training data from {train_csv}...")
    try:
        df = pd.read_csv(train_csv)
    except Exception as e:
        print(f"Error loading {train_csv}: {e}")
        return

    if 'label' not in df.columns or 'smiles' not in df.columns:
        print("Error: Dataset must contain 'smiles' and 'label' columns.")
        # Fallback for flexibility
        if 'delta_g' in df.columns: df['label'] = df['delta_g'] 
        if 'final_deltag' in df.columns: df['label'] = df['final_deltag']

    # Drop missing
    df = df.dropna(subset=['smiles', 'label'])
    
    # 2. Generate Features (SMILES -> ECFP4)
    X, valid_idx = generate_fingerprints(df['smiles'].values)
    y = df['label'].values[valid_idx]
    
    print(f"Valid Samples: {len(X)} / {len(df)}")
    
    # 3. Split Data
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
    
    # 4. Train Model
    print("Training XGBoost Regressor...")
    model = xgb.XGBRegressor(
        objective='reg:squarederror',
        n_estimators=500,
        learning_rate=0.05,
        max_depth=6,
        n_jobs=-1
    )
    model.fit(X_train, y_train)
    
    # 5. Evaluate
    preds_train = model.predict(X_train)
    preds_test = model.predict(X_test)
    
    mse_train = mean_squared_error(y_train, preds_train)
    mse_test = mean_squared_error(y_test, preds_test)
    r2_test = r2_score(y_test, preds_test)
    
    print(f"Train RMSE: {np.sqrt(mse_train):.4f}")
    print(f"Test RMSE:  {np.sqrt(mse_test):.4f}")
    print(f"Test R^2:   {r2_test:.4f}")
    
    # 6. Save Model
    print(f"Saving model to {output_model}...")
    os.makedirs(os.path.dirname(output_model), exist_ok=True)
    joblib.dump(model, output_model)
    print("Model saved successfully.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--train", required=True, help="Path to Training CSV (smiles, label)")
    parser.add_argument("--output_model", required=True, help="Path to save model (e.g. model.joblib)")
    # Removed prediction args for training-only job simplicity
    
    args = parser.parse_args()
    
    train_xgboost(args.train, args.output_model)

