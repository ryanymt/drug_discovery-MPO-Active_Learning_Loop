#!/bin/bash

# Optimized Gnina Docking Script (PARALLEL VERSION)
# Changes:
# 1. Splits input SMILES into 4 chunks.
# 2. Runs 4 concurrent Gnina pipelines (obabel processing + gnina scoring) in background.
# 3. Merges results.
# 4. Periodic sync from the main merger (conceptually hard with splits) or just sync final.
#    Actually, we can have each worker write to its own CSV, and we sync the directory?
#    Or just wait for all to finish. Since this is "Cleanup", speed is key.

INPUT_SMILES_FILE=$1
INPUT_PROTEIN_FILE=$2
OUTPUT_DIR=$3
TASK_INDEX=$4
EXHAUSTIVENESS=${5:-1} # Default low for screening, can raise for cleanup if needed

WORK_DIR="/tmp/gnina_work_${TASK_INDEX}"
mkdir -p "${WORK_DIR}"
cd "${WORK_DIR}"

# 1. Prepare Protein (Once)
cp "${INPUT_PROTEIN_FILE}" ./protein.pdb
if ! obabel protein.pdb -O protein.pdbqt -xr; then
    echo "CRITICAL ERROR: Failed to convert protein to PDBQT"
    exit 1
fi

PROTEIN_PDBQT="protein.pdbqt"

# 2. Split Input
# Detect Available CPUs
CORES=$(nproc)
echo "Detected ${CORES} CPU cores. Launching ${CORES} parallel workers."

# Split into N chunks (numeric suffixes)
split -n l/${CORES} -d "${INPUT_SMILES_FILE}" chunk_

# 3. Define Worker Function
run_worker() {
    local chunk_file=$1
    local worker_id=$2
    local output_csv="results_${worker_id}.csv"
    local log_file="worker_${worker_id}.log"
    
    echo "worker_${worker_id} started on ${chunk_file}" > ${log_file}
    
    # Process loop (Simulating the original script's per-mol logic)
    # We can reuse the logic, but streamlined.
    
    i=0
    while IFS= read -r smiles || [[ -n "$smiles" ]]; do
        if [ -z "$smiles" ]; then continue; fi
        
        # Unique IDs per worker
        LIGAND_PDBQT="ligand_${worker_id}_${i}.pdbqt"
        DOCKED_PDBQT="docked_${worker_id}_${i}.pdbqt"
        
        # 1. Obabel
        if ! obabel -:"${smiles}" -O ${LIGAND_PDBQT} --gen3d -p 7.4 --partialcharge gasteiger >> ${log_file} 2>&1; then
            echo "${smiles},NaN,Obabel_Failed" >> ${output_csv}
            i=$((i+1))
            continue
        fi
        
        if [ ! -s "${LIGAND_PDBQT}" ]; then
            echo "${smiles},NaN,PDBQT_Empty" >> ${output_csv}
            i=$((i+1))
            continue
        fi
        
        # 2. Gnina
        # Note: All workers share the GPU. L4 can handle context switching, but we might saturate VRAM?
        # L4 has 24GB. Tiny molecules take almost nothing. Should be fine.
        if ! gnina --receptor ${PROTEIN_PDBQT} \
              --ligand ${LIGAND_PDBQT} \
              --autobox_ligand ${LIGAND_PDBQT} \
              --autobox_add 8 \
              --out ${DOCKED_PDBQT} \
              --cpu 1 \
              --exhaustiveness ${EXHAUSTIVENESS} >> ${log_file} 2>&1; then
             # gnina failed
             true
        fi
        
        # 3. Score
        SCORE=$(grep -A 1 "CNN pose score" ${log_file} | tail -n 1 | awk '{print $4}')
        if [ -z "$SCORE" ]; then
            SCORE="NaN"
            STATUS="Gnina_Failed"
        else
            STATUS="Success"
        fi
        
        echo "${smiles},${SCORE},${STATUS}" >> ${output_csv}
        i=$((i+1))
        
    done < "${chunk_file}"
    
    echo "worker_${worker_id} finished" >> ${log_file}
}

# 4. Launch Workers
pids=""
for f in chunk_*; do
    # Extract suffix as ID (00, 01, 02...)
    id=${f#chunk_}
    echo "Launching Worker $id on $f"
    run_worker "$f" "$id" &
    pids="$pids $!"
done

# 5. Wait
for pid in $pids; do
    wait $pid
done

# 6. Merge
FINAL_OUTPUT_CSV="${OUTPUT_DIR}/predictions_${TASK_INDEX}.csv"
FINAL_LOG_FILE="${OUTPUT_DIR}/docking_log_${TASK_INDEX}.txt"

echo "smiles,docking_score,status" > ${FINAL_OUTPUT_CSV}
cat results_*.csv >> ${FINAL_OUTPUT_CSV}
cat worker_*.log > ${FINAL_LOG_FILE}

echo "Parallel Docking Complete."
