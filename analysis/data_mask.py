import netCDF4 as nc4
import numpy as np
from scipy.spatial import cKDTree
import xarray as xr

# def get_ds(chlor_file_path, var = 'chla'):
#     ds = nc4.Dataset(chlor_file_path)
    
#     print(ds)
#     print(ds.dims)
#     print(ds.data_vars)
#     # get the mean through time
#     # ds = ds.mean(dim='time')
#     x = np.array(ds.variables['x'][:]).astype(float)
#     y = np.array(ds.variables['y'][:]).astype(float)
#     # data = ds.variables[var][:, :]
#     # data = ds through time
#     data = ds.mean(dim='time').variables[var][:, :]

#     data = np.transpose(data)
#     return x, y, data

def get_ds(chlor_file_path, var):
    ds = xr.open_dataset(chlor_file_path)
    x = ds['x'].values.astype(float)
    y = ds['y'].values.astype(float)
    data = ds[var].mean(dim='time')
    data = data.transpose()
    return x, y, data

def make_distance_mask(data, x, y):
    # Disko Bay bounds (e.g.)
    min_x = -463350
    min_y = -2364239
    max_x = -115044
    max_y = -2049628

    # Read in the BedMachine mask
    ds = nc4.Dataset('/Users/tara/Documents/SJSU/MLML/Grids/BedMachine/BedMachineGreenland-v5.nc')
    bm_x = ds.variables['x'][:]
    bm_y = ds.variables['y'][:]
    bm_mask = ds.variables['mask'][:, :]
    ds.close()

    # Subset to Disko Bay
    x_indices = np.logical_and(bm_x>min_x, bm_x<max_x)
    y_indices = np.logical_and(bm_y>min_y, bm_y<max_y)
    bm_x = np.array(bm_x[x_indices]).astype(float)
    bm_y = np.array(bm_y[y_indices]).astype(float)
    bm_mask = bm_mask[y_indices, :]
    bm_mask = bm_mask[:, x_indices]

    # interpolate BedMachine to the grid as a starting point
    bm_X, bm_Y = np.meshgrid(bm_x, bm_y)
    X, Y = np.meshgrid(x, y)

    # find the points corresponding to ice/land
    land_indices = bm_mask.ravel()>0
    land_x = bm_X.ravel()[land_indices]
    land_y = bm_Y.ravel()[land_indices]

    # make a mask to fill in
    mask = np.zeros_like(data)

    # Create a KDTree from the land points
    tree = cKDTree(np.column_stack((land_x, land_y)))

    # Calculate the distance to the closest land point for each point in the grid
    distances, _ = tree.query(np.column_stack((X.ravel(), Y.ravel())))

    # Reshape the distances array to the shape of the mask
    mask = distances.reshape(mask.shape)

    return mask

def mask_data(data, mask, threshold=1000):
    # define a distance threshold
    # adjust this as needed to eliminate most of the noise
    distance_threshold = threshold # 2km = 2000m

    # define the mask
    distance_mask = mask<distance_threshold

    # mask the data
    data_masked = np.copy(data)
    data_masked[distance_mask] = np.nan

    return data_masked