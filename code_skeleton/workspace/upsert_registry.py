
import pandas as pd
from google.cloud import bigquery
import argparse
import sys
import os

def upsert_registry(csv_path, project_id, dataset_id, table_id="molecule_registry"):
    """
    Upserts molecules from CSV into BigQuery Registry using MERGE.
    1. Load CSV to Staging Table.
    2. MERGE Staging -> Registry.
    3. Drop Staging.
    """
    client = bigquery.Client(project=project_id)
    
    # 1. Config
    staging_table_id = f"{table_id}_staging"
    dataset_ref = client.dataset(dataset_id)
    registry_ref = dataset_ref.table(table_id)
    staging_ref = dataset_ref.table(staging_table_id)
    
    print(f"--- Upserting {csv_path} to {project_id}.{dataset_id}.{table_id} ---")
    
    # 2. Load to Staging
    job_config = bigquery.LoadJobConfig(
        source_format=bigquery.SourceFormat.CSV,
        skip_leading_rows=1,
        autodetect=True,
        write_disposition=bigquery.WriteDisposition.WRITE_TRUNCATE # Replace staging
    )

    with open(csv_path, "rb") as source_file:
        job = client.load_table_from_file(source_file, staging_ref, job_config=job_config)

    job.result()  # Wait for job to complete
    print(f"Loaded {job.output_rows} rows to staging table.")

    # 3. Perform MERGE
    # Matches on 'molecule_hash'
    # When NOT MATCHED -> INSERT
    # When MATCHED -> UPDATE (e.g. update run_id or just ignore? Let's UPDATE to keep latest metadata if needed, or ignore. User asked for unique. Let's Ignore duplicates.)
    
    merge_query = f"""
    MERGE `{project_id}.{dataset_id}.{table_id}` T
    USING `{project_id}.{dataset_id}.{staging_table_id}` S
    ON T.molecule_hash = S.molecule_hash
    WHEN NOT MATCHED THEN
      INSERT (molecule_hash, smiles, global_id, sdf_gcs_path, run_id, created_at)
      VALUES (S.molecule_hash, S.smiles, S.global_id, S.sdf_gcs_path, S.run_id, S.created_at)
    """
    # Note: IF ignoring updates, we don't add WHEN MATCHED clause.
    # This ensures "First Writer Wins" or "Existing Stays". 
    # If we want "Latest Writer Wins", we'd add update.
    # Given the user complained about duplicates, preserving the existing unique ID is safest.

    query_job = client.query(merge_query)
    query_job.result()
    print("MERGE operation complete.")
    
    # 4. Cleanup
    client.delete_table(staging_ref, not_found_ok=True)
    print("Staging table dropped.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv", required=True)
    parser.add_argument("--project", default="gcda-apac-sc")
    parser.add_argument("--dataset", default="bioops_platform")
    args = parser.parse_args()
    
    if not os.path.exists(args.csv):
        print(f"File not found: {args.csv}")
        sys.exit(1)
        
    upsert_registry(args.csv, args.project, args.dataset)
