import os
import glob
import xarray as xr

def get_files(seawifs_folder, modis_folder, cell):
    seawifs_files = os.path.join(seawifs_folder, cell)
    modis_files = os.path.join(modis_folder, cell)

    seawifs_files = glob.glob(os.path.join(seawifs_files, '*.nc'))
    modis_files = glob.glob(os.path.join(modis_files, '*.nc'))

    seawifs_files.sort()
    modis_files.sort()

    return seawifs_files, modis_files

def coarsen_data(files, coarsen_factor):
    coarse_data = xr.open_mfdataset(files)
    coarse_data = coarse_data.coarsen(x=coarsen_factor, y=coarsen_factor).mean()
    return coarse_data

def combine_arrays(mean_seawifs, mean_modis):
    combined_array = xr.concat([mean_seawifs, mean_modis], dim='time')
    return combined_array

def write_file(save_path, cell, combined_array):
    save_file = os.path.join(save_path, cell+'_combined.nc')
    if not os.path.isfile(save_file):
        combined_array.to_netcdf(save_file)
    else:
        print('WARNING: File already exists for cell: ', cell)
    return