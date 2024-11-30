import os
import pandas as pd
import glob

# Find all learning_log.csv files
# Experiment runs downloaded from https://zenodo.org/records/13984543
file_paths = glob.glob(
    os.path.expanduser(
        "~/Downloads/experiment_runs/grid_search_compound/*/learning_log.csv"
    )
)
# List to hold DataFrames
dataframes = []

for file_path in file_paths:
    # Extract the first 8 characters of the folder name
    folder_name = os.path.basename(os.path.dirname(file_path))

    # Read the CSV file
    df = pd.read_csv(file_path)

    # Add a new column for the folder name
    df["folder"] = folder_name

    # Append to the list
    dataframes.append(df)

# Concatenate all DataFrames into one
combined_df = pd.concat(dataframes, ignore_index=True)

# Save the combined DataFrame to a new CSV file
combined_df.to_csv("output/leave_compound_out_learning_log.csv.gz", index=False)
