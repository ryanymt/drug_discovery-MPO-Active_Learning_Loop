#!/bin/bash
set -e

# --- CONFIG ---
PROJECT_ID="gcda-apac-sc"
DATASET="bioops_platform"
LOCATION="us-central1"

echo "Using Project: $PROJECT_ID"
echo "Dataset: $DATASET"

# 1. Create Dataset
if ! bq ls --project_id $PROJECT_ID | grep -q $DATASET; then
    echo "Creating dataset..."
    bq --project_id $PROJECT_ID mk --dataset --location=$LOCATION $DATASET
else
    echo "Dataset $DATASET already exists."
fi

# 2. Create Tables (using Schema Definitions inline or from JSON)

# Table: molecule_registry
echo "Creating/Updating molecule_registry..."
bq mk --table \
  --schema "molecule_hash:STRING,smiles:STRING,global_id:STRING,sdf_gcs_path:STRING,run_id:STRING,created_at:TIMESTAMP" \
  --time_partitioning_field created_at \
  $PROJECT_ID:$DATASET.molecule_registry || true

# Table: screening_results
echo "Creating/Updating screening_results..."
bq mk --table \
  --schema "molecule_hash:STRING,docking_score:FLOAT,qed_score:FLOAT,sa_score:FLOAT,logp:FLOAT,mw:FLOAT,toxicity_label:INTEGER,final_score:FLOAT,score_method:STRING" \
  --clustering_fields molecule_hash \
  $PROJECT_ID:$DATASET.screening_results || true

# Table: final_affinity
echo "Creating/Updating final_affinity..."
bq mk --table \
  --schema "molecule_hash:STRING,affinity_score:FLOAT,method:STRING,uncertainty:FLOAT,stage:STRING" \
  $PROJECT_ID:$DATASET.final_affinity || true

# View: v_leaderboard
echo "Creating View: v_leaderboard..."
QUERY="
CREATE OR REPLACE VIEW \`$PROJECT_ID.$DATASET.v_leaderboard\` AS
SELECT
    r.smiles,
    a.affinity_score,
    a.method,
    s.docking_score,
    s.qed_score,
    a.molecule_hash
FROM \`$PROJECT_ID.$DATASET.molecule_registry\` r
JOIN \`$PROJECT_ID.$DATASET.final_affinity\` a ON r.molecule_hash = a.molecule_hash
LEFT JOIN \`$PROJECT_ID.$DATASET.screening_results\` s ON r.molecule_hash = s.molecule_hash
ORDER BY a.affinity_score ASC
"
bq --project_id $PROJECT_ID query --use_legacy_sql=false "$QUERY"

echo "BigQuery Setup Complete."
