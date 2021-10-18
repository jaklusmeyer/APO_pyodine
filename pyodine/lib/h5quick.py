# A few simple tools for quick reading of HDF5-data

import h5py


def h5print(h, level=0):
    """Recursively print the structure of a HDF5 file or group
    
    Args:
        h (h5py file or group): The h5py object.
        level (Optional[int]): The level within the structure at which to start
            (needed to print the whole structure recursively).
    """
    if level == 0:
        print(h.name)
    for k in h.keys():
        print('    ' * level + '└ ' + k)
        if isinstance(h[k], h5py.Group):
            h5print(h[k], level + 1)


def h5data(h):
    """Retrieve a HDF5 dataset or group
    
    If the handle is a group, return value is a dict of datasets (recursive).
    
    Args:
        h (h5py dataset or group): The h5py to return.
    
    Return:
        dict or dataset: The data packed into a dictionary resembling the
            structure of the h5py object.
    """
    if isinstance(h, h5py.Dataset):
        return h[()]
    else:
        return {k: h5data(h[k]) for k in h.keys()}


def h5get(filename, item):
    """Retrieve a named item from a HDF5 file
    
    If the item is a group, return value is a dict of datasets (recursive).
    
    Args:
        filename (str): The path to the HDF5 file.
        item (str): A key to the dataset or group of interest.
    
    Return:
        dict or dataset: The data, either as dictionary or a dataset.
    """
    with h5py.File(filename, 'r') as h:
        return h5data(h[item])


# Helper function
def dict_to_group(my_dict, base_group, new_group_name):
    """Pack data into a h5py object
    
    Create a hdf5 group with name `new_group_name` in the existing
    group handle `base_group` (could be the root) and fill in
    named datasets from `my_dict`
    
    Args:
        my_dict (dict): The data to pack.
        base_group (h5py group): The group handle to pack the data into.
        new_group_name (str): The name of the group.
    """
    group = base_group.create_group(new_group_name)
    for k in my_dict:
        group[k] = my_dict[k]