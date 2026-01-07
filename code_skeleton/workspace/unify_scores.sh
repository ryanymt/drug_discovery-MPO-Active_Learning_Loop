#!/bin/bash
set -e

PROJECT_ID="gcda-apac-sc"
DATASET_ID="bioops_platform"
PREDICTIONS_CSV="gs://ryanymt/output/inference/predictions.csv"

SCREENING_TABLE_SQL="$PROJECT_ID.$DATASET_ID.screening_results"
STAGING_TABLE_SQL="$PROJECT_ID.$DATASET_ID.predictions_staging"
GROMACS_TABLE_SQL="$PROJECT_ID.$DATASET_ID.final_affinity"

# CLI Format (project:dataset.table)
STAGING_TABLE_CLI="$PROJECT_ID:$DATASET_ID.predictions_staging"

echo "--- Unifying Scores in $SCREENING_TABLE_SQL ---"

# 1. Load XGBoost Predictions into Staging
echo "Loading predictions into staging table..."
bq load --source_format=CSV --skip_leading_rows=1 --autodetect --replace \
  $STAGING_TABLE_CLI $PREDICTIONS_CSV

# 2. Update with XGBoost Scores (Baseline)
echo "Merging XGBoost Scores..."
QUERY_XGB="
UPDATE \`$SCREENING_TABLE_SQL\` T
SET 
    final_score = S.final_score,
    score_method = 'xgboost'
FROM \`$STAGING_TABLE_SQL\` S
WHERE T.molecule_hash = S.molecule_hash
"
bq query --use_legacy_sql=false "$QUERY_XGB"

# 3. Validation: Compare Predictions vs Ground Truth (The 1k)
echo "--- Model Validation (Pred vs Actual) ---"
QUERY_VAL="
SELECT 
    CORR(S.final_score, G.delta_g) as correlation,
    SQRT(AVG(POW(S.final_score - G.delta_g, 2))) as rmse,
    COUNT(*) as count
FROM \`$STAGING_TABLE_SQL\` S
JOIN \`$GROMACS_TABLE_SQL\` G ON S.molecule_hash = G.molecule_hash
WHERE NOT IS_NAN(G.delta_g)
"
bq query --use_legacy_sql=false "$QUERY_VAL"

# 4. Overwrite with GROMACS Scores (Ground Truth)
echo "Merging GROMACS Scores (Overwriting Proxy)..."
QUERY_GMX="
UPDATE \`$SCREENING_TABLE_SQL\` T
SET 
    final_score = G.delta_g,
    score_method = 'gromacs'
FROM \`$GROMACS_TABLE_SQL\` G
WHERE T.molecule_hash = G.molecule_hash
  AND G.delta_g IS NOT NULL
  AND NOT IS_NAN(G.delta_g)
"
bq query --use_legacy_sql=false "$QUERY_GMX"

# 5. Clean up
echo "Cleaning up staging..."
bq rm -f -t $STAGING_TABLE_CLI


echo "Done."
