#!/bin/bash
# set -e  <-- DISABLED for robustness

# --- Defaults ---
MODE="all"
INPUT_DIR="/mnt/disks/gcs/input"
OUTPUT_DIR="/mnt/disks/gcs/output"
MDP_DIR="/opt/fep"
PROTEIN_PDB="${INPUT_DIR}/protein.pdb"
LIGAND_PDB="${INPUT_DIR}/ligand.pdb"

# --- Argument Parsing ---
while [[ $# -gt 0 ]]; do
    case $1 in
        --run) MODE="run" ;;
        --id) JOB_ID=$2; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

if [ -z "$JOB_ID" ]; then
    JOB_ID="Unknown"
fi

source /opt/gromacs/bin/GMXRC

# --- Setup Logging ---
WORK_DIR=$(mktemp -d)
cd ${WORK_DIR}

# Log everything to a specific file that we sync
DEBUG_LOG="${WORK_DIR}/debug_${JOB_ID}.log"
RESULTS_FILE="${OUTPUT_DIR}/${JOB_ID}_results.csv"  # GCS Destination
LOCAL_RESULTS_FILE="${WORK_DIR}/results.csv"        # Local Buffer
# Sync Log File Destination
GCS_LOG_FILE="${OUTPUT_DIR}/debug_${JOB_ID}.log"

# Function to log and sync
log_msg() {
    echo "$(date) - $1" | tee -a $DEBUG_LOG
    cp $DEBUG_LOG $GCS_LOG_FILE
}

log_msg "--- STARTING PROBE JOB ${JOB_ID} ---"
log_msg "Checking Hardware:"
nvidia-smi >> $DEBUG_LOG 2>&1 || log_msg "nvidia-smi failed"
gmx --version >> $DEBUG_LOG 2>&1 || log_msg "gmx failed"

python3 -c "import rdkit; print('RDKit Available')" || (log_msg "Error: RDKit not found" && exit 1)

# Check Inputs
# Check Inputs
BATCH_CSV="${INPUT_DIR}/${JOB_ID}.csv"
SDF_INPUT="${INPUT_DIR}/ligand.sdf"

if [ -f "$BATCH_CSV" ]; then
    log_msg "Found Batch CSV: $BATCH_CSV"
elif [ -f "$SDF_INPUT" ]; then
    log_msg "Found Batch SDF: $SDF_INPUT. Proceeding without CSV."
    # Create valid dummy variable to prevent python script error if it references it
    BATCH_CSV="/dev/null"
else
    log_msg "Error: Neither Batch CSV ($BATCH_CSV) nor SDF ($SDF_INPUT) found."
    ls -lR $INPUT_DIR >> $DEBUG_LOG
    exit 1
fi

# Determine MDP Source (with logging)
if [ -f "${INPUT_DIR}/minim_stable.mdp" ]; then
    MDP_SOURCE="${INPUT_DIR}"
    log_msg "Using Custom MDPs from ${MDP_SOURCE}"
else
    MDP_SOURCE="${MDP_DIR}"
    log_msg "WARNING: Custom MDPs not found. Using Container Defaults from ${MDP_SOURCE}"
fi

# DEBUG: Print nsteps
log_msg "Checking nsteps in production mdp:"
grep "nsteps" ${MDP_SOURCE}/prod_stable.mdp >> $DEBUG_LOG 2>&1 || log_msg "Could not grep prod_stable.mdp"

# Init Results
echo "SMILES,DeltaG,Status" > $LOCAL_RESULTS_FILE
cp $LOCAL_RESULTS_FILE $RESULTS_FILE
cp $DEBUG_LOG $GCS_LOG_FILE

# Generate helper scripts
cat <<EOF > generate_pdbs.py
import csv
import sys
import os
from rdkit import Chem
from rdkit.Chem import AllChem

csv_file = "${BATCH_CSV}"
os.makedirs("mols", exist_ok=True)
valid_mols = []

# --- NEW: SDF INPUT SUPPORT ---
sdf_input = "${INPUT_DIR}/ligand.sdf"
if os.path.exists(sdf_input):
    # Found matching SDF, utilizing valid coordinates.
    suppl = Chem.SDMolSupplier(sdf_input)
    for i, mol in enumerate(suppl):
        if mol is None: continue
        mol_id = f"mol_{i}"
        # SDF has 3D coords, so we DO NOT embed or optimize
        mol = Chem.AddHs(mol, addCoords=True)
        valid_mols.append({"id": mol_id, "mol": mol})
else:
    # Fallback to CSV + Embedding
    if os.path.exists(csv_file):
        with open(csv_file, 'r') as f:
            reader = csv.DictReader(f)
            for i, row in enumerate(reader):
                smiles = row.get('SMILES') or row.get('smiles')
                if not smiles: continue
                mol_id = f"mol_{i}"
                
                try:
                    mol = Chem.MolFromSmiles(smiles)
                    if not mol: continue
                    mol = Chem.AddHs(mol)
                    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
                    AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
                    
                    valid_mols.append({"id": mol_id, "mol": mol})
                except:
                    pass

for m in valid_mols:
    mol_id = m["id"]
    out_path = f"mols/{mol_id}/ligand.pdb"
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    Chem.MolToPDBFile(m["mol"], out_path)
    print(f"{mol_id} {out_path}")
    
    # FORCE RENAMING OF RESIDUE 'UNL' TO 'MOL'
    # This ensures consistency through Acpype -> GROMACS
    os.system(f"sed -i 's/UNL/MOL/g' {out_path}")
    os.system(f"sed -i 's/LIG/MOL/g' {out_path}") # Safety catch for other defaults
EOF

python3 generate_pdbs.py > mol_list.txt

# Loop
while read mol_id pdb_path; do
    log_msg ">>> STARTING PIPELINE for $mol_id <<<"
    
    MOL_DIR="${WORK_DIR}/run_${mol_id}"
    mkdir -p $MOL_DIR
    cd $MOL_DIR
    
    cp $WORK_DIR/$pdb_path ./ligand.pdb
    cp $PROTEIN_PDB ./protein.pdb
    
    # 1. Setup
    if ! python3 /usr/local/bin/fep_setup.py --protein protein.pdb --ligand ligand.pdb >> $DEBUG_LOG 2>&1; then
        log_msg "Failed Setup for $mol_id"
        echo "$mol_id,NaN,Failed_Setup" >> $LOCAL_RESULTS_FILE
        cp $LOCAL_RESULTS_FILE $RESULTS_FILE
        cp $DEBUG_LOG $GCS_LOG_FILE
        continue
    fi
    
    # 2. Minim
    if ! gmx grompp -f ${MDP_SOURCE}/minim_stable.mdp -c system.gro -p topol.top -o minim.tpr -maxwarn 2 >> $DEBUG_LOG 2>&1 || \
       ! gmx mdrun -v -deffnm minim -ntmpi 1 >> $DEBUG_LOG 2>&1; then
        log_msg "Failed Minim for $mol_id"
        echo "$mol_id,NaN,Failed_Minim" >> $LOCAL_RESULTS_FILE
        cp $LOCAL_RESULTS_FILE $RESULTS_FILE
        cp $DEBUG_LOG $GCS_LOG_FILE
        continue
    fi
    
    # 3. Equil
    if ! gmx grompp -f ${MDP_SOURCE}/nvt_stable.mdp -c minim.gro -r minim.gro -p topol.top -o nvt.tpr -maxwarn 2 >> $DEBUG_LOG 2>&1 || \
       ! gmx mdrun -v -deffnm nvt -ntmpi 1 -nb gpu >> $DEBUG_LOG 2>&1; then
        log_msg "Failed Equil for $mol_id"
        echo "$mol_id,NaN,Failed_Equil" >> $LOCAL_RESULTS_FILE
        cp $LOCAL_RESULTS_FILE $RESULTS_FILE
        cp $DEBUG_LOG $GCS_LOG_FILE
        continue
    fi
    
    log_msg "Starting Production (GPU Check)..."
    
    # 4. Prod
    if ! gmx grompp -f ${MDP_SOURCE}/prod_stable.mdp -c nvt.gro -p topol.top -o production.tpr -maxwarn 2 >> $DEBUG_LOG 2>&1 || \
       ! gmx mdrun -v -deffnm production -ntmpi 1 -nb gpu >> $DEBUG_LOG 2>&1; then
        log_msg "Failed Prod for $mol_id"
        echo "$mol_id,NaN,Failed_Production" >> $LOCAL_RESULTS_FILE
        cp $LOCAL_RESULTS_FILE $RESULTS_FILE
        cp $DEBUG_LOG $GCS_LOG_FILE
        continue
    fi
    
    log_msg "Starting Analysis..."
    
    # 5. Analysis (FIXED: PBC & FITTING)
    
    # --- PHASE 1: Center Protein & Unwrap (-pbc mol) ---
    # Fixes the "Vanishing Ligand" / "DeltaG=0.00" issue.
    log_msg "S5: Centering Protein (Protein System)..."
    echo "Protein System" | gmx trjconv -s production.tpr -f production.xtc -o production_centered.xtc -center -pbc mol >> $DEBUG_LOG 2>&1
    
    # --- PHASE 2: Fit to Reference (-fit rot+trans) ---
    # Aligns tumbling complex for stable analysis.
    log_msg "S5: Fitting Trajectory (Backbone System)..."
    echo "Backbone System" | gmx trjconv -s production.tpr -f production_centered.xtc -o production_fit.xtc -fit rot+trans >> $DEBUG_LOG 2>&1
    
    # Safety Check
    if [ ! -f "production_fit.xtc" ]; then
        log_msg "CRITICAL ERROR: production_fit.xtc was not created."
        echo "$mol_id,NaN,Failed_PBC_Fix" >> $LOCAL_RESULTS_FILE
        continue
    fi
    
    # (Continue with standard mmpbsa.in setup below...)
    cat <<IN_EOF > mmpbsa_dynamic.in
    
    
&general
sys_name="ProtLig",
startframe=1,
endframe=9999,
/
&gb
igb=5, saltcon=0.150,
/
# &nmode
# nmstartframe=1, nmendframe=9999,
# nmode_istrng=0.150,
# /
IN_EOF
    
    # Create Index File
    # Ensure MOL is recognized as a group.
    # The residue name in fep_setup.py is "MOL".
    # gmx make_ndx auto-selects residues into groups.
    echo -e "q\n" | gmx make_ndx -f production.tpr -o index.ndx >> $DEBUG_LOG 2>&1
    
    # Run MMPBSA with Index File
    # Use -ci for index file (MMPBSA specific flag) and -cg for group selection.
    # We rely on topol.top for ligand parameters (no -lm needed).
    echo "--- Running gmx_MMPBSA for $mol_id ---"
    
    # Check for MPI and parallelize
    # CRITICAL: mpi4py IS REQUIRED for gmx_MMPBSA to run with mpirun.
    # If mpi4py is missing, mpirun will launch N independent serial jobs that clobber each other (Deadlock).
    HAS_MPI4PY=false
    if python3 -c "from mpi4py import MPI" &> /dev/null; then
        HAS_MPI4PY=true
    else
        log_msg "mpi4py not found. Attempting to install..."
        if python3 -m pip install mpi4py; then
            HAS_MPI4PY=true
            log_msg "mpi4py installed successfully."
        else
            log_msg "Failed to install mpi4py. Will fallback to serial or multiprocessing if available."
        fi
    fi
    # SUPPRESS GUI and fix environment
    export QT_QPA_PLATFORM=offscreen
    export GMX_MMPBSA_NO_GUI=1
    
    # Fix for QStandardPaths: XDG_RUNTIME_DIR not set (Safety catch)
    export XDG_RUNTIME_DIR=/tmp/runtime-root
    mkdir -p $XDG_RUNTIME_DIR
    chmod 700 $XDG_RUNTIME_DIR

    # AUTO-KILL GUI WATCHDOG
    # gmx_MMPBSA_ana hangs in headless environments even without -O flag.
    # We launch a background sniper to kill it if it appears.
    (
        while true; do
            if pgrep -f "gmx_MMPBSA_ana" > /dev/null; then
                log_msg "WATCHDOG: GUI Analyzer detected. Killing it to prevent hang..."
                pkill -f "gmx_MMPBSA_ana"
            fi
            sleep 5
        done
    ) &
    WATCHDOG_PID=$!

    if [ "$HAS_MPI4PY" = true ] && command -v mpirun &> /dev/null; then
        NCORES=$(nproc)
        log_msg "MPI/mpi4py Detected. Running on $NCORES cores."
        mpirun --allow-run-as-root --use-hwthread-cpus -np $NCORES gmx_MMPBSA -i mmpbsa_dynamic.in -cs production.tpr -ct production_fit.xtc -cp topol.top -ci index.ndx -cg Protein MOL >> $DEBUG_LOG 2>&1 || true
    else
        log_msg "MPI or mpi4py Not Found. Running serially."
        gmx_MMPBSA -i mmpbsa_dynamic.in -cs production.tpr -ct production_fit.xtc -cp topol.top -ci index.ndx -cg Protein MOL >> $DEBUG_LOG 2>&1 || true
    fi
    
    # Kill the watchdog
    kill $WATCHDOG_PID 2>/dev/null
    
    # DEBUG: List generated files to confirm output name
    log_msg "--- Files in $MOL_DIR ---"
    ls -l >> $DEBUG_LOG
    
    # PARSING STRATEGY:
    # 1. Try FINAL_RESULTS_MMPBSA.dat (cleanest)
    # 2. Try _GMXMMPBSA_info (raw python dict/text)
    # 3. Try parsing the log footer (last resort)
    
    DG_VAL=""
    
    # Attempt 1: Standard File
    if [ -f "FINAL_RESULTS_MMPBSA.dat" ]; then
         cat FINAL_RESULTS_MMPBSA.dat >> $DEBUG_LOG
         # Look for "Delta G binding = -30.50 +/- 3.25"
         DG_VAL=$(grep "Delta G binding" FINAL_RESULTS_MMPBSA.dat | awk -F'=' '{print $2}' | awk '{print $1}')
    fi
    
    # Attempt 2: Raw Stats (if DG_VAL is empty or 0.00 which is suspicious)
    if [ -z "$DG_VAL" ] || [ "$DG_VAL" == "0.00" ]; then
        if [ -f "_GMXMMPBSA_info" ]; then
            log_msg "Attempting to parse _GMXMMPBSA_info..."
             DG_VAL=$(grep "GENERAL" _GMXMMPBSA_info | grep "TOTAL" | tail -n 1 | awk '{print $NF}')
        fi
    fi
    
    # Attempt 3: Log Sniffing (Reliable if printed to stdout)
    # FIX: Exclude "Attempting" to avoid matching the log message itself
    if [ -z "$DG_VAL" ] || [ "$DG_VAL" == "0.00" ]; then
         log_msg "Attempting to sniff log table for TOTAL..."
         # Extract the 'ΔTOTAL' line from the log we just wrote/captured
         # Line format in log: "ΔTOTAL                   -29.50          3.25"
         # Use grep -v to ignore our own log messages (which start with Date or contains "Attempting")
         DG_VAL=$(grep "ΔTOTAL" $DEBUG_LOG | grep -v "Attempting" | grep -v "Sat Jan" | tail -n 1 | awk '{print $2}')
    fi

    if [ -z "$DG_VAL" ]; then
        DG_VAL="NaN"
        STATUS="Failed_Analysis"
        log_msg "Analysis Failed: DeltaG not found or parse error"
    else
        STATUS="Success"
        log_msg "SUCCESS: $mol_id DeltaG = $DG_VAL"
    fi
    
    echo "$mol_id,$DG_VAL,$STATUS" >> $LOCAL_RESULTS_FILE
    
    # Sync
    cp $LOCAL_RESULTS_FILE $RESULTS_FILE
    cp $DEBUG_LOG $GCS_LOG_FILE
    
    rm -f *.trr *.edr minim* nvt*
    cd ${WORK_DIR}

done < mol_list.txt

log_msg "Rank Verification Complete"
cp $DEBUG_LOG $GCS_LOG_FILE
sync
