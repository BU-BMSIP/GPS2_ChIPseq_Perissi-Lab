import pandas as pd

# Read sample data
samples = pd.read_csv("Sample_Data.csv")

# Ensure unique columns
if samples.columns.duplicated().any():
    raise ValueError("Duplicate columns found in Sample_Data.csv")

# Convert DataFrame to dictionary using 'File name' as the key
sample_dict = samples.set_index('File name').T.to_dict()

# Save preprocessed sample dictionary to a JSON file for easy access
pd.DataFrame(sample_dict).to_json("samples.json")