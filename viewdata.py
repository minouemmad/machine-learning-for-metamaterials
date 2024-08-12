import h5py

import numpy as np

# Disable truncation of arrays in the output
np.set_printoptions(threshold=np.inf)
filename = "C:\\Users\\maemmad\\Desktop\\machine-learning-for-metamaterials-1\\data_generation\\Data\\data_rte_gen2lay5mat_100n_v-tma.h5"

def print_structure(name, obj):
    """Helper function to recursively print the structure of the HDF5 file."""
    print(f"{name}: {type(obj)}")
    if isinstance(obj, h5py.Group):
        for key, value in obj.attrs.items():
            print(f"  Attribute: {key} = {value}")
    elif isinstance(obj, h5py.Dataset):
        print(f"  Shape: {obj.shape}")
        for key, value in obj.attrs.items():
            print(f"  Attribute: {key} = {value}")

with h5py.File(filename, "r") as f:
    # Print the structure of the file
    f.visititems(print_structure)
    
    # Get the first object name/key
    a_group_key = list(f.keys())[0]
    
    # Load the dataset
    ds_obj = f[a_group_key]
    ds_arr = ds_obj[()]  # Load the data into a numpy array
    
    # Print a summary of the data
    print(f"\nDataset shape: {ds_arr.shape}")
    print(f"Dataset type: {ds_arr.dtype}")
    print(f"Dataset type: {list(f.keys())[0]}")
    
    # Check for and print dataset attributes (titles, labels, etc.)
    print("\nDataset Attributes:")
    for key, value in ds_obj.attrs.items():
        print(f"  {key}: {value}")
    
    # Print the first few rows to get a sense of the data
    # print("\nFirst 5 rows of the dataset:")
    # print(ds_arr[2])
    # print(len(ds_arr[2]))

    # Print titles if they exist
    if 'titles' in ds_obj.attrs:
        print("\nTitles (if available):")
        print(ds_obj.attrs['titles'])

