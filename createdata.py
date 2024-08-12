import h5py
import numpy as np

# Define the data you want to store
data = np.array([
    [25., 45., 65., 4.74504213, 4.74378326, 4.74252941],
    [25., 45., 65., 5.3463119, 5.34272373, 5.33914025],
    [25., 45., 65., 4.35533835, 4.35253365, 4.34974097],
    [25., 45., 65., 4.61989077, 4.61669694, 4.61351522],
    [25., 45., 65., 5.63853309, 5.63422717, 5.62993195]
])

# Create a new HDF5 file
filename = "C:\\Users\\maemmad\\Desktop\\testdata.h5"
with h5py.File(filename, "w") as f:
    # Create a dataset within the file
    dset = f.create_dataset("data", data=data)
    
    # Optionally add attributes to the dataset (e.g., titles, units)
    dset.attrs['description'] = "Sample dataset with angles and reflectance data"
    dset.attrs['angles'] = "90 degrees"
    dset.attrs['measurement'] = "Reflectance data"

print(f"HDF5 file '{filename}' created successfully.")
