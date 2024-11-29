import json
import glob
import pandas as pd
from pathlib import Path
import os

# Option 1: Using os.path.expanduser
folders = glob.glob(
    os.path.expanduser("~/Downloads/experiment_runs/grid_search_compound/*/")
)

results = []
for folder in folders:
    try:
        # Read all three files
        with open(Path(folder) / "config.json") as f:
            config = json.load(f)
        with open(Path(folder) / "test_scores.json") as f:
            test = json.load(f)
        with open(Path(folder) / "val_scores.json") as f:
            val = json.load(f)

        # Combine into single dict
        row = {**config}  # Start with config
        row.update({f"test_{k}": v for k, v in test.items()})  # Add test scores
        row.update({f"val_{k}": v for k, v in val.items()})  # Add val scores
        results.append(row)

    except Exception as e:
        print(f"Error processing {folder}: {e}")

# Convert to DataFrame and save
pd.DataFrame(results).to_csv("output/leave_compound_out_combined.csv", index=False)
print(f"Processed {len(results)} folders")
