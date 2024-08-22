import h5py
import pandas as pd
import os

def extract_and_save(dataset, path=""):
    # If the dataset is a group, we need to explore it further
    if isinstance(dataset, h5py.Group):
        for name, item in dataset.items():
            extract_and_save(item, f"{path}/{name}")
    elif isinstance(dataset, h5py.Dataset):
        # If it's a dataset, convert it to a DataFrame and save it as a CSV file
        
        # Ensure the directory exists
        os.makedirs(os.path.dirname(path), exist_ok=True)
        
        # Convert the dataset to a DataFrame
        df = pd.DataFrame(dataset[()])
        
        # Save the DataFrame to a CSV file
        df.to_csv(f'{path}.csv', index=False)
    else:
        print(f"Unknown type encountered at {path}")

# Open the HDF5 file
with h5py.File('mse_space\\dielectric_nk_data.h5', 'r') as h5file:
    for name, dataset in h5file.items():
        extract_and_save(dataset, name)
