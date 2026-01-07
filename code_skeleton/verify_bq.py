from google.cloud import bigquery
import os

project_id = "lifescience-project-469915"
dataset_id = "durg_discovery_results"
table_id = "docking_admet_results"

client = bigquery.Client(project=project_id)

dataset_ref = client.dataset(dataset_id)

print(f"Checking dataset: {dataset_id}...")
try:
    client.get_dataset(dataset_ref)
    print("Dataset exists.")
except Exception as e:
    print(f"Dataset not found or access denied: {e}")
    exit(1)

print(f"Checking table: {table_id}...")
try:
    tables = list(client.list_tables(dataset_ref))
    found = False
    for table in tables:
        print(f"Found table: {table.table_id}")
        if table.table_id == table_id:
            found = True
    
    if found:
        print(f"SUCCESS: Table {table_id} exists.")
        
        # Fetch a few rows
        table_ref = dataset_ref.table(table_id)
        rows = client.list_rows(table_ref, max_results=5)
        print("\nSample Data:")
        for row in rows:
            print(dict(row))
            
    else:
        print(f"FAILURE: Table {table_id} does NOT exist in dataset.")

except Exception as e:
    print(f"Error listing tables: {e}")
