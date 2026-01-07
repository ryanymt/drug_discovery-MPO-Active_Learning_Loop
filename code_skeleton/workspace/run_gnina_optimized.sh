#!/bin/bash

# Optimized Gnina Docking Script (Robust + Periodic Sync)
# Changes:
# 1. Dynamic CPU usage
# 2. Configurable Exhaustiveness
# 3. ROBUST Error Handling (No set -e)
# 4. PERIODIC SYNC to GCS (Every 10 molecules)

# --- GPU Environment Setup ---
export LD_LIBRARY_PATH=/usr/local/nvidia/lib64:/usr/local/cuda/lib64:${LD_LIBRARY_PATH}
nvidia-smi

# --- Args ---
INPUT_SMILES_FILE=$1
INPUT_PROTEIN_FILE=$2
OUTPUT_DIR=$3
TASK_INDEX=$4
EXHAUSTIVENESS=${5:-1}

# --- Setup ---
WORK_DIR="/tmp/gnina_work_${TASK_INDEX}"
mkdir -p "${WORK_DIR}"
cd "${WORK_DIR}"

# --- File Paths ---
PROTEIN_PDBQT="protein.pdbqt"
RESULTS_CSV="results.csv"
FINAL_OUTPUT_CSV="${OUTPUT_DIR}/predictions_${TASK_INDEX}.csv"
FINAL_LOG_FILE="${OUTPUT_DIR}/docking_log_${TASK_INDEX}.txt"

echo "Starting ROBUST OPTIMIZED GNINA docking task ${TASK_INDEX}..."
echo "Configuration: Exhaustiveness=${EXHAUSTIVENESS}, CPUs=$(nproc)"

# Prepare protein
cp "${INPUT_PROTEIN_FILE}" ./protein.pdb
if ! obabel protein.pdb -O ${PROTEIN_PDBQT} -xr; then
    echo "CRITICAL ERROR: Failed to convert protein to PDBQT"
    exit 1
fi

if [ ! -s "${PROTEIN_PDBQT}" ]; then
    echo "CRITICAL ERROR: Protein PDBQT empty"
    exit 1
fi

echo "smiles,docking_score,status" > ${RESULTS_CSV}

# Process each SMILES string
i=0
while IFS= read -r smiles || [[ -n "$smiles" ]]; do
    if [ -z "$smiles" ]; then
        continue
    fi

    LIGAND_PDBQT="ligand_${i}.pdbqt"
    DOCKED_PDBQT="docked_${i}.pdbqt"
    LOG_FILE="log_${i}.txt"

    echo "Processing SMILES ${i}: ${smiles}" > ${LOG_FILE}

    # 1. Convert SMILES to 3D PDBQT
    if ! obabel -:"${smiles}" -O ${LIGAND_PDBQT} --gen3d -p 7.4 --partialcharge gasteiger >> ${LOG_FILE} 2>&1; then
        echo "Error: Obabel failed" >> ${LOG_FILE}
        echo "${smiles},NaN,Obabel_Failed" >> ${RESULTS_CSV}
        i=$((i+1))
        continue
    fi

    if [ ! -s "${LIGAND_PDBQT}" ]; then
        echo "Error: PDBQT empty" >> ${LOG_FILE}
        echo "${smiles},NaN,PDBQT_Empty" >> ${RESULTS_CSV}
        i=$((i+1))
        continue
    fi

    # 2. Run GNINA docking
    if ! gnina --receptor ${PROTEIN_PDBQT} \
          --ligand ${LIGAND_PDBQT} \
          --autobox_ligand ${LIGAND_PDBQT} \
          --autobox_add 8 \
          --out ${DOCKED_PDBQT} \
          --log ${LOG_FILE} \
          --cpu $(nproc) \
          --exhaustiveness ${EXHAUSTIVENESS} >> ${LOG_FILE} 2>&1; then
          
          echo "Warning: Gnina returned non-zero exit code" >> ${LOG_FILE}
    fi

    # 3. Extract Score
    SCORE=$(grep -A 1 "CNN pose score" ${LOG_FILE} | tail -n 1 | awk '{print $4}')
    
    if [ -z "$SCORE" ]; then
        SCORE="NaN"
        STATUS="Gnina_Failed"
    else
        STATUS="Success"
    fi

    # 4. Append
    echo "${smiles},${SCORE},${STATUS}" >> ${RESULTS_CSV}

    # 5. Periodic Sync (Every 10 molecules) -- ADDED
    if (( i % 10 == 0 )); then
        echo "Syncing results to GCS (Index $i)..."
        cp ${RESULTS_CSV} ${FINAL_OUTPUT_CSV}
    fi

    i=$((i+1))
done < "${INPUT_SMILES_FILE}"

# --- Finalization ---
echo "Docking complete. Copying results."
mkdir -p "${OUTPUT_DIR}"
cp ${RESULTS_CSV} ${FINAL_OUTPUT_CSV}
find . -name "log_*.txt" -exec cat {} + > ${FINAL_LOG_FILE}
