import json
from kfp import dsl
from kfp.compiler import Compiler
from google.cloud import aiplatform

# Import the component
from batch_component import submit_batch_job

# --- Project Configuration ---
PROJECT_ID = "gcda-apac-sc"
REGION = "us-central1"
PIPELINE_ROOT = "gs://ryanymt/output/pipeline-root-v4"
BASE_GCS_PATH = "/mnt/disks/gcs"

# --- Spec Loading Helper ---
def load_and_modify_spec(filename, modifier_func=None, **kwargs):
    with open(filename, "r") as f:
        spec = json.load(f)
    if modifier_func:
        modifier_func(spec, **kwargs)
    return spec

# --- Modifiers ---
def modify_pocket2mol(spec):
    container = spec['taskGroups'][0]['taskSpec']['runnables'][0]['container']
    container['imageUri'] = "us-central1-docker.pkg.dev/gcda-apac-sc/drug-discovery-containers/pocket2mol:v3"
    cmds = container['commands']
    cmd_str = cmds[2]
    
    # Configure for V4 loop
    # 500 samples for quick iteration
    cmd_str = cmd_str.replace("num_samples: 100000", "num_samples: 500") 
    cmd_str = cmd_str.replace("num_samples: 1000", "num_samples: 500")

    # Output directory
    new_outdir = f"{BASE_GCS_PATH}/output/v4_loop_1"
    
    # Environment
    conda_setup = "source /opt/conda/etc/profile.d/conda.sh && conda activate Pocket2Mol && "
    
    # Args
    import re
    cmd_str = re.sub(r"--outdir\s+\S+", f"--outdir {new_outdir}", cmd_str)
    # Use internal 4yhj reference
    cmd_str = re.sub(r"--pdb_path\s+\S+", " --pdb_path ./example/4yhj.pdb", cmd_str)
    # Set threshold
    if "--rejection_threshold" not in cmd_str:
        cmd_str += " --rejection_threshold -5.0"
    
    # Ensure Conda
    if "conda activate" not in cmd_str:
        cmd_str = conda_setup + cmd_str
        
    cmds[2] = cmd_str
    
def modify_filters(spec, input_flag):
    # Point to the Pocket2Mol output locations
    # Pocket2Mol outputs SMILES.txt to the outdir
    shared_input_path = f"{BASE_GCS_PATH}/input/v4_loop_1/SMILES.txt"
    # Note: We rely on the fact that /mnt/disks/gcs/input and .../output map to the same bucket
    # So writing to output/v4_loop_1 makes it available at input/v4_loop_1 in the next step
    # IF the next step mounts 'ryanymt/output' to /mnt/disks/gcs/input
    
    cmds = spec['taskGroups'][0]['taskSpec']['runnables'][0]['container']['commands']
    
    if input_flag: 
        # Replace flag logic
        for i, token in enumerate(cmds):
            if input_flag in token:
                 import re
                 cmds[i] = re.sub(fr"{input_flag}\s+\S+", f"{input_flag} {shared_input_path}", token)
    else:
        # Gnina logic (positional)
        cmd_str = cmds[2]
        parts = cmd_str.split()
        if "/usr/local/bin/run_gnina_docking.sh" in parts:
             idx = parts.index("/usr/local/bin/run_gnina_docking.sh")
             if idx + 1 < len(parts):
                 parts[idx+1] = shared_input_path
             cmds[2] = " ".join(parts)

def modify_joiner(spec):
    # Joiner needs paths to all previous outputs
    # Inputs:
    smiles_in = f"{BASE_GCS_PATH}/input/v4_loop_1/SMILES.txt"
    rdkit_in = f"{BASE_GCS_PATH}/input/v4_loop_1/rdkit_scores.csv" # RDKit usually outputs to cwd or specified
    # We need to ensure Filters output to v4_loop_1 too!
    # ... For this MVP, let's assume Filters are configured to output to standard names in the input dir or similar
    # Actually, the task-specs for filters might default to specific paths.
    # Let's override them in modify_filters? 
    # For now, let's assume standard names in the same folder for simplicity of this orchestration script
    
    # Just setting up the joiner command args
    container = spec['taskGroups'][0]['taskSpec']['runnables'][0]['container']
    # The joiner container likely has an ENTRYPOINT or CMD we need to overwrite or append to.
    # Assuming the joiner image has `python join_results.py` as default or we invoke it.
    
    cmd = [
        "python3", "join_results.py",
        "--smiles", smiles_in,
        "--rdkit", f"{BASE_GCS_PATH}/input/v4_loop_1/rdkit_scores.csv",
        "--txgemma", f"{BASE_GCS_PATH}/input/v4_loop_1/txgemma_results.csv",
        "--gnina", f"{BASE_GCS_PATH}/input/v4_loop_1/gnina_out", # Directory
        "--output", f"{BASE_GCS_PATH}/output/v4_loop_1/joined_results.csv",
        "--batch_id", "v4_loop_1"
    ]
    container['commands'] = cmd
    # Ensure image
    container['imageUri'] = "us-central1-docker.pkg.dev/gcda-apac-sc/drug-discovery-containers/joiner:v1"

def modify_selection(spec):
    # Run Vizier
    container = spec['taskGroups'][0]['taskSpec']['runnables'][0]['container']
    cmd = [
        "python3", "run_vizier_optimization.py",
        "--output", f"{BASE_GCS_PATH}/output/v4_loop_1/selected_candidates.csv"
    ]
    container['commands'] = cmd
    # Using the same image as joiner or a general python image? 
    # Joiner image has pandas/google-cloud, should be enough.
    container['imageUri'] = "us-central1-docker.pkg.dev/gcda-apac-sc/drug-discovery-containers/joiner:v1"

def modify_fep_setup(spec):
    container = spec['taskGroups'][0]['taskSpec']['runnables'][0]['container']
    # Use selected candidates
    # setup script needs to know where to find them. 
    # If run_fep.sh --setup uses a hardcoded input, we might need to change it or move files.
    # Let's assume we pass the candidates file or it looks in a standard place.
    # For now, we'll copy the selected file to 'ligands.csv' or similar if needed, 
    # but `run_fep.sh` in setup mode likely processes a specific list.
    # We'll just run it as is, assuming it picks up from the standard input dir where we wrote 'selected_candidates.csv'
    
    container['commands'] = ["/bin/bash", "-c", "/usr/local/bin/run_fep.sh --setup"]
    container['imageUri'] = "us-central1-docker.pkg.dev/gcda-apac-sc/drug-discovery-containers/gromacs-fep:v1"

def modify_fep_run(spec):
    container = spec['taskGroups'][0]['taskSpec']['runnables'][0]['container']
    container['commands'] = ["/bin/bash", "-c", "/usr/local/bin/run_fep.sh --run $BATCH_TASK_INDEX"]
    container['imageUri'] = "us-central1-docker.pkg.dev/gcda-apac-sc/drug-discovery-containers/gromacs-fep:v1"
    
    spec['taskGroups'][0]['taskCount'] = "13"
    spec['taskGroups'][0]['parallelism'] = "13"
    
    spec['allocationPolicy']['instances'][0]['policy']['machineType'] = "g2-standard-4"
    spec['allocationPolicy']['instances'][0]['installGpuDrivers'] = True
    spec['allocationPolicy']['instances'][0]['policy']['accelerators'] = [
        {"type": "nvidia-l4", "count": "1", "installGpuDrivers": True}
    ]

def modify_fep_analyze(spec):
    container = spec['taskGroups'][0]['taskSpec']['runnables'][0]['container']
    container['commands'] = ["/bin/bash", "-c", "/usr/local/bin/run_fep.sh --analyze"]
    container['imageUri'] = "us-central1-docker.pkg.dev/gcda-apac-sc/drug-discovery-containers/gromacs-fep:v1"

# --- Pipeline Definition ---
@dsl.pipeline(
    name="drug-discovery-v4-active-learning",
    description="V4: Generation -> Filtering -> Join -> Selection -> FEP"
)
def discovery_pipeline_v4():
    
    # 1. Generation
    p2m_spec = load_and_modify_spec("pocket2mol-task-spec.json", modify_pocket2mol)
    p2m = submit_batch_job(
        project=PROJECT_ID, location=REGION, job_spec=json.dumps(p2m_spec)
    ).set_display_name("1. Generation (Pocket2Mol)")
    
    # 2. Filtering
    gnina_spec = load_and_modify_spec("gnina-task-spec.json", modify_filters, input_flag=None)
    gnina = submit_batch_job(
        project=PROJECT_ID, location=REGION, job_spec=json.dumps(gnina_spec)
    ).set_display_name("2a. Filter (Gnina)").after(p2m)

    tx_spec = load_and_modify_spec("txgemma-task-spec.json", modify_filters, input_flag="--input-file")
    txgemma = submit_batch_job(
        project=PROJECT_ID, location=REGION, job_spec=json.dumps(tx_spec)
    ).set_display_name("2b. Filter (TxGemma)").after(p2m)

    rdkit_spec = load_and_modify_spec("rdkit-task-spec.json", modify_filters, input_flag="--input_csv")
    rdkit = submit_batch_job(
        project=PROJECT_ID, location=REGION, job_spec=json.dumps(rdkit_spec)
    ).set_display_name("2c. Filter (RDKit)").after(p2m)
    
    # 3. Joiner
    join_spec = load_and_modify_spec("join-task-spec.json", modify_joiner)
    joiner = submit_batch_job(
        project=PROJECT_ID, location=REGION, job_spec=json.dumps(join_spec)
    ).set_display_name("3. Data Joiner").after(gnina, txgemma, rdkit)
    
    # 4. Selection
    # Re-use join-task-spec structure as base for python script execution
    sel_spec = load_and_modify_spec("join-task-spec.json", modify_selection)
    selection = submit_batch_job(
        project=PROJECT_ID, location=REGION, job_spec=json.dumps(sel_spec)
    ).set_display_name("4. Selection (Vizier)").after(joiner)
    
    # 5. FEP Workflow
    # A. Setup
    fep_setup_spec = load_and_modify_spec("gromacs-task-spec.json", modify_fep_setup)
    fep_setup = submit_batch_job(
        project=PROJECT_ID, location=REGION, job_spec=json.dumps(fep_setup_spec)
    ).set_display_name("5a. FEP Setup").after(selection)
    
    # B. Run
    fep_run_spec = load_and_modify_spec("gromacs-task-spec.json", modify_fep_run)
    fep_run = submit_batch_job(
        project=PROJECT_ID, location=REGION, job_spec=json.dumps(fep_run_spec)
    ).set_display_name("5b. FEP Run (Parallel)").after(fep_setup)
    
    # C. Analyze
    fep_analyze_spec = load_and_modify_spec("gromacs-task-spec.json", modify_fep_analyze)
    fep_analyze = submit_batch_job(
        project=PROJECT_ID, location=REGION, job_spec=json.dumps(fep_analyze_spec)
    ).set_display_name("5c. FEP Analysis").after(fep_run)

if __name__ == "__main__":
    Compiler().compile(pipeline_func=discovery_pipeline_v4, package_path="drug_discovery_pipeline_v4.json")
    
    aiplatform.init(project=PROJECT_ID, location=REGION)
    job = aiplatform.PipelineJob(
        display_name="v4-active-learning-loop",
        template_path="drug_discovery_pipeline_v4.json",
        pipeline_root=PIPELINE_ROOT,
        enable_caching=False
    )
    # job.submit() # Uncomment to run
    print("V4 Pipeline Compiled Successfully.")
