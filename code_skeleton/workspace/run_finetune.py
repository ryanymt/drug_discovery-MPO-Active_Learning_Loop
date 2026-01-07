import os
import sys
import subprocess
import shutil
import tarfile
import random
import torch
import pickle
from google.cloud import storage

def download_blob(bucket_name, source_blob_name, destination_file_name):
    """Downloads a blob from the bucket."""
    try:
        storage_client = storage.Client()
        bucket = storage_client.bucket(bucket_name)
        blob = bucket.blob(source_blob_name)
        blob.download_to_filename(destination_file_name)
        print(f"Blob {source_blob_name} downloaded to {destination_file_name}.")
    except Exception as e:
        print(f"Error downloading {source_blob_name}: {e}")
        raise e

def patch_dataset_loader():
    """
    Downloads patched pl.py and overwrites the container's version.
    """
    print("Patching dataset loader for multiprocessing...")
    try:
        # Download patched file
        download_blob("ryanymt", "input/scripts/pl_patched.py", "pl_patched.py")
        
        # Find where the library is installed
        # We need to be careful here. The original code might be in site-packages
        # or relative to the current working directory if it's a local import.
        # Let's try to locate "utils.datasets.pl" by importing it from current dir structure if possible,
        # or rely on the fact that we are in the root.
        
        # Actually, we can just search for the file if we are in the root source dir.
        # But let's rely on finding the file via the python path.
        # Since we can't easily import "utils" without conflict if we are not careful...
        # We will try a simpler approach script-side: 
        # find the file "pocket2mol-retrain/utils/datasets/pl.py" or similar in standard paths.
        
        # Simpler: The training script runs 'import utils.datasets'. 
        # So we should be able to find 'utils/datasets/pl.py' relative to CWD if we are in the repo root?
        # The Docker container WORKDIR is usually the repo root.
        # Let's try to overwrite "utils/datasets/pl.py" AND "pocket2mol_retrain/utils/datasets/pl.py" if they exist.
        
        possible_paths = [
            "utils/datasets/pl.py",
            "pocket2mol-retrain/utils/datasets/pl.py",
            "/app/pocket2mol-retrain/utils/datasets/pl.py",
            "/opt/conda/envs/Pocket2Mol/lib/python3.8/site-packages/pocket2mol_retrain/utils/datasets/pl.py"
        ]
        
        patched = False
        for path in possible_paths:
            if os.path.exists(path):
                print(f"Found pl.py at {path}, overwriting...")
                shutil.copy("pl_patched.py", path)
                patched = True
                
        # Also try via import as a fallback
        if not patched:
            try:
                import importlib.util
                # Try likely module names
                for mod_name in ['utils.datasets.pl', 'pocket2mol_retrain.utils.datasets.pl']:
                    try:
                        spec = importlib.util.find_spec(mod_name)
                        if spec and spec.origin:
                            print(f"Found via module {mod_name}: {spec.origin}")
                            shutil.copy("pl_patched.py", spec.origin)
                            patched = True
                            break # Found and patched, no need to check other modules
                    except Exception:
                        pass # Ignore import errors for non-existent modules
            except Exception as e:
                print(f"Error trying to resolve via import: {e}")
                
        if patched:
            print("Successfully patched pl.py for multiprocessing.")
        else:
            print("WARNING: Could not find pl.py to patch. Using slower single-threaded loader.")
            
    except Exception as e:
        print(f"WARNING: Failed to patch dataset loader: {e}")
        print("Using default single-threaded loader.")

def extract_tar(tar_path, dest_dir):
    print(f"Extracting {tar_path} to {dest_dir}...")
    with tarfile.open(tar_path, 'r:gz') as tar:
        tar.extractall(path=dest_dir)
    print("Extraction complete.")

def create_split_file(data_dir, split_path):
    print(f"Generating split file at {split_path}...")
    files = [f for f in os.listdir(data_dir) if f.endswith('.sdf')]
    print(f"Found {len(files)} SDF files.")
    
    random.shuffle(files)
    split_idx = int(len(files) * 0.9) # 90% train, 10% val
    
    train_files = files[:split_idx]
    val_files = files[split_idx:]
    
    # Remove extension for P2M dataset index usually? 
    # Pocket2Mol usually expects list of (pdb, ligand) tuples or just names if 'pl' dataset.
    # If the dataloader loads by name from the directory.
    # It usually expects the name WITHOUT extension if it constructs paths like {name}/{name}_ligand.sdf
    # But our directory is flat: {hash}.sdf.
    # This might be an issue if the dataset class expects a specific folder structure (e.g. crossdocked format).
    # Crossdocked format: {pair_id}/{pair_id}_pocket.pdb, {pair_id}/{pair_id}_ligand.sdf
    # We only have SDFs (ligands). 
    # Does P2M need the pocket? YES. It generates IN A POCKET.
    # CRITICAL: We forgot the pockets!
    # The elite candidates come from generated molecules, which were generated in *some* pocket.
    # Usually the SAME pocket for all (if relying on single-target generation).
    # Task 2264 (01_Pocket2Mol_Generator.md): "Input: Protein Pocket Structure (.pdb)".
    # If we are fine-tuning on generated molecules (Experience Replay), we must provide the protein context they were generated for.
    # In Phase 1/2, we used `4yhj.pdb` (or `AF-P32298.pdb`).
    # If all 100k were generated for `4yhj.pdb`, then we just need to provide that single PDB for all.
    # The dataset loader might expect {name}_pocket.pdb for each entry.
    # Solution: We will duplicate the reference PDB for each training sample, OR modify the dataset code (hard).
    # Duplicating 10k PDBs is wasteful but reliable if we can't change code.
    # Or symlink?
    
    print("Mocking pocket PDBs...")
    ref_pdb = "input/protein.pdb" # We will ensure this exists
    if not os.path.exists(ref_pdb):
        print("WARNING: Reference PDB not found at input/protein.pdb!")
    
    for f in files:
        name = os.path.splitext(f)[0]
        # P2M Dataset (PL) often expects:
        # root/name/name_ligand.sdf
        # root/name/name_pocket.pdb
        # We need to restructure our flat directory to this if using standard loader.
        # Or if it supports flat: root/name_ligand.sdf, root/name_pocket.pdb.
        
        # Let's assume standard Crossdocked split logic:
        # split = {'train': ['pair1', 'pair2'], ...}
        
    split_dict = {'train': train_files, 'valid': val_files, 'test': []}
    torch.save(split_dict, split_path)
    print("Split file saved.")

def fix_environment():
    # Reuse v7 hot fix logic
    print("--- HOTFIX: Force Re-installing Compatible PyTorch Geometric ---")
    try:
        subprocess.check_call([sys.executable, "-m", "pip", "uninstall", "-y", "torch-scatter", "torch-sparse", "torch-cluster", "torch-spline-conv", "torch-geometric"])
        subprocess.check_call([
            sys.executable, "-m", "pip", "install",
            "torch-scatter==2.0.9", "torch-sparse==0.6.12", "torch-cluster==1.5.9", "torch-spline-conv==1.2.1",
            "-f", "https://data.pyg.org/whl/torch-1.10.0+cu113.html"
        ])
        subprocess.check_call([sys.executable, "-m", "pip", "install", "torch-geometric==2.0.4", "more-itertools"])
    except Exception as e:
        print(f"Hotfix failed: {e}")

def run_finetune():
    print("--- Starting Pocket2Mol Fine-Tuning (V8: Elite Data) ---")
    
    fix_environment()
    
    # Setup Dirs
    os.makedirs("configs", exist_ok=True)
    os.makedirs("ckpt", exist_ok=True)
    os.makedirs("data/finetune", exist_ok=True)
    os.makedirs("logs", exist_ok=True)
    os.makedirs("input", exist_ok=True) # For PDB
    
    # Download Inputs
    # 1. Config
    download_blob("ryanymt", "input/scripts/finetune.yml", "configs/finetune.yml")
    # 2. Base Model
    download_blob("ryanymt", "input/models/pretrained.pt", "ckpt/pocket2mol.pt")
    # 3. Training Data (Tarball)
    download_blob("ryanymt", "input/training/elite_sdfs/elite_sdfs.tar.gz", "elite_sdfs.tar.gz")
    # 4. Reference Protein (CRITICAL)
    download_blob("ryanymt", "input/AF-P32298.pdb", "input/protein.pdb") # Assuming this is the target
    
    # Extract
    extract_tar("elite_sdfs.tar.gz", "data/finetune")
    
    # Restructure Data (Recursive Discovery)
    print("Structuring Processed Data...")
    
    # Use os.walk to find SDFs anywhere in data/finetune
    sdf_files = []
    for root, dirs, files in os.walk("data/finetune"):
        for f in files:
            if f.endswith(".sdf") and not f.endswith("_ligand.sdf"):
                sdf_files.append(os.path.join(root, f))
    
    print(f"Found {len(sdf_files)} SDF files (recursive scan).")
    
    names = []
    ref_pdb_path = os.path.abspath("input/protein.pdb")
    
    # We need to flatten them effectively for the split file?
    # Or just ensure they are named correctly in place.
    # P2M likely expects a flat directory if we provide `path: ./data/finetune`.
    # If files are in `data/finetune/elite_prep/x.sdf`, we should MOVE them to `data/finetune`.
    
    for src in sdf_files:
        filename = os.path.basename(src)
        name = os.path.splitext(filename)[0]
        names.append(name)
        
        # Target paths in the ROOT of data/finetune
        dst_sdf = os.path.join("data/finetune", f"{name}_ligand.sdf")
        dst_pdb = os.path.join("data/finetune", f"{name}_pocket.pdb")
        
        # Move/Rename SDF
        if src != dst_sdf:
            # Check if we are overwriting itself (e.g. flat dir)
            # If flat: src = data/finetune/idx.sdf. dst = data/finetune/idx_ligand.sdf. Safe.
            # If nested: src = data/finetune/sub/idx.sdf. dst = data/finetune/idx_ligand.sdf. Safe.
            shutil.move(src, dst_sdf)
            
        # Symlink Protein
        if not os.path.exists(dst_pdb):
            os.symlink(ref_pdb_path, dst_pdb)
            
    # Cleanup empty dirs if nested
    # (Optional)
    
    print(f"Prepared {len(names)} training samples.") 
                
    # Create Split
    # P2M split file contains list of NAMES (without suffixes)
    
    # 0. Patch the dataset loader (Hotfix for parallel loading)
    patch_dataset_loader()

    # 1. Environment Fix (Required for PyG / CUDA)
    # List of (pocket_fn, ligand_fn, _, rmsd_str)
    index = []
    split_names = [] # List of tuples for split file
    
    for name in names:
        pocket_fn = f"{name}_pocket.pdb"
        ligand_fn = f"{name}_ligand.sdf"
        index.append((pocket_fn, ligand_fn, 0, 0))
        split_names.append((pocket_fn, ligand_fn))
        
    with open("data/finetune/index.pkl", "wb") as f:
        pickle.dump(index, f)
    print(f"Generated data/finetune/index.pkl with {len(index)} samples.")

    # Create Split (using tuples to match name2id)
    # The keys in dataset.name2id are (protein_filename, ligand_filename)
    split_dict = {
        'train': split_names[:int(len(split_names)*0.9)],
        'valid': split_names[int(len(split_names)*0.9):],
        'test': []
    }
    torch.save(split_dict, "data/finetune/split.pt")
    
    # Config is ready at configs/finetune.yml (downloaded)
    print("Using config: configs/finetune.yml")

    # Run Training
    cmd = [
        sys.executable, "/app/train_finetune.py",
        "--config", "configs/finetune.yml",
        "--resume", "ckpt/pocket2mol.pt",
        "--logdir", "logs"
    ]
    
    # Handle LD_LIBRARY_PATH (from v7)
    print("--- Fixing LD_LIBRARY_PATH ---")
    lib_path = find_library("libcusparse.so.11")
    
    env = os.environ.copy()
    
    if lib_path:
        lib_dir = os.path.dirname(lib_path)
        print(f"Found libcusparse.so.11 at: {lib_path}")
        
        current_ld = env.get("LD_LIBRARY_PATH", "")
        if current_ld:
             new_ld = f"{current_ld}:{lib_dir}"
        else:
             new_ld = lib_dir
             
        env["LD_LIBRARY_PATH"] = new_ld
        print(f"Updated LD_LIBRARY_PATH (Appended): {new_ld}")
    else:
        print("WARNING: Could not find libcusparse.so.11! Training will likely fail.")

    # Executing Training
    print("--- Executing Training ---")
    subprocess.check_call(cmd, env=env)
    
    # Upload Logs
    upload_folder("ryanymt", "logs", "output/active_learning/cycle2_model/logs")
    
    # Upload Checkpoints (CRITICAL FIX)
    upload_folder("ryanymt", "ckpt", "output/active_learning/cycle2_model/ckpt")
    
    # Showcase Vertex AI Metadata
    log_vertex_metadata(
        dataset_uri="gs://ryanymt/input/training/elite_sdfs/elite_sdfs.tar.gz",
        model_uri="gs://ryanymt/output/active_learning/cycle2_model/pocket2mol.pt", # Approximate final loc
        params={"epochs": 10, "batch_size": 4, "lr": 1e-4}
    )

def find_library(name):
    print(f"Searching for {name} in /opt/conda...")
    for root, dirs, files in os.walk("/opt/conda"):
        if name in files:
            return os.path.join(root, name)
    return None

def log_vertex_metadata(dataset_uri, model_uri, params):
    print("--- Logging to Vertex AI Metadata ---")
    try:
        # Install if missing
        subprocess.check_call([sys.executable, "-m", "pip", "install", "google-cloud-aiplatform", "--quiet"])
        from google.cloud import aiplatform
        
        aiplatform.init(project="gcda-apac-sc", location="us-central1")
        
        # Start an Execution
        with aiplatform.start_execution(
            schema_title="system.ContainerExecution", 
            display_name="pocket2mol-finetune-v8"
        ) as execution:
            
            # Log Parameters
            print(f"Logging params: {params}")
            # execution.log_params(params) # Not direct method on execution context in all SDK versions, usually attached to Run. 
            # For Metadata directly:
            # We create Artifacts.
            
            # Input Artifact
            ds_artifact = aiplatform.Artifact.create(
                schema_title="system.Dataset",
                display_name="elite_sdfs_cycle_2",
                uri=dataset_uri
            )
            execution.assign_input_artifacts([ds_artifact])
            
            # Output Artifact
            model_artifact = aiplatform.Artifact.create(
                schema_title="system.Model",
                display_name="pocket2mol_cycle_2_model",
                uri=model_uri
            )
            execution.assign_output_artifacts([model_artifact])
            
            print("Vertex AI Metadata logged successfully.")
            
    except Exception as e:
        print(f"Warning: Failed to log Vertex Metadata: {e}")

def upload_folder(bucket_name, source_folder, prefix):
    storage_client = storage.Client()
    bucket = storage_client.bucket(bucket_name)
    for root, dirs, files in os.walk(source_folder):
        for file in files:
            local_path = os.path.join(root, file)
            rel_path = os.path.relpath(local_path, source_folder)
            blob_path = os.path.join(prefix, rel_path)
            bucket.blob(blob_path).upload_from_filename(local_path)

if __name__ == "__main__":
    run_finetune()
