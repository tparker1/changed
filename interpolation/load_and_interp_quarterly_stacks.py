import netCDF4 as nc
import toolbox.grid_generation as grid
import toolbox.initiation.grid_construction as grid_construction
import numpy as np
import toolbox.reprojection as reprojection
from scipy.interpolate import griddata
import xarray as xr
import os
import glob
import dask


def get_lat_lon(file):
    """
    Get the lat and lon from a file
    :param file: file path to the netcdf file
    :return: lons, lats
    """
    ds = nc.Dataset(file)
    lons = ds.variables['lon'][:]
    lats = ds.variables['lat'][:]
    ds.close()
    return lons, lats


def load_file(file):
    """
    Load a netcdf file
    :param file: file path to the netcdf file
    :return: chla, time
    """
    ds = nc.Dataset(file)
    chla = ds.variables['chlor_a'][:]
    time = ds.variables['time'][:]
    ds.close()

    return chla, time


def load_file_xr(file):
    """
    Load a netcdf file using xarray
    :param file: file path to the netcdf file
    :return: chla, time
    """
    ds = xr.open_dataset(file)
    chla = ds['chlor_a'].values
    time = ds['time'].values
    ds.close()

    return chla, time


# def specify_grid(CC, grid_of_interest):
#     """
#     Specify the grid of interest. ie. the grid you want interpolation for.
#     :param CC: Chla Object
#     :param grid_of_interest: grid of interest as a list: [min_row, min_col, max_row, max_col]
#     """
#     grid.store_cells(CC, grid_of_interest)
#     return


def buffer_extents_and_create_grid(CC, buffer_multiplier=100):
    """
    Buffer the extents of the grid by 100 cells and create a grid from the buffered extents.
    :param CC: Chla Object
    :param buffer_multiplier: how much buffer to add. Default = 100. (buffer = CC.posting * buffer_multiplier)
    """
    buffer = CC.posting * buffer_multiplier

    CC.buffered_extents = [
        CC.extents[0] - buffer, CC.extents[1] - buffer, CC.extents[2] + buffer,
        CC.extents[3] + buffer
    ]
    grid_construction.create_grid_from_extents(CC)
    return


def reproject(lons, lats):
    """
    Reproject data to the desired projection
    :param lons: lons from the netcdf file
    :param lats: lats from the netcdf file
    """
    inputCRS = 4326
    outputCRS = 3413

    xx, yy = np.meshgrid(lons, lats)

    polygon_array = np.array([xx.flatten(), yy.flatten()]).T

    output_polygon = reprojection.reproject_polygon(polygon_array, inputCRS,
                                                    outputCRS)

    x = output_polygon[:, 0].reshape(xx.shape)
    y = output_polygon[:, 1].reshape(yy.shape)

    return x, y


def interpolate_data(CC, x, y, chla):
    """
    Interpolate the data to the grid of interest
    :param CC: Chla Object
    :param x: x coordinates in 3413
    :param y: y coordinates in 3413
    :param chla: chla data
    :return: interpolated_chla
    """
    # Assuming you have your data loaded as x, y, and chla
    # Also assuming you have chl_grid_x and chl_grid_y defined

    # Define the bounds of your cell in 3413]
    minx, miny, maxx, maxy = CC.buffered_extents

    # Crop the data to the bounds of your cell
    mask = (x >= minx) & (x <= maxx) & (y >= miny) & (y <= maxy)
    x_cropped = x[mask]
    y_cropped = y[mask]
    chla_cropped = chla[:, mask]

    # Create a meshgrid for the cell
    xx, yy = np.meshgrid(CC.chl_grid_x, CC.chl_grid_y)

    # Reshape the meshgrid to match the dimensions of the data
    xx = xx.flatten()
    yy = yy.flatten()

    # Initialize an empty array for the interpolated data
    interpolated_chla = np.zeros(
        (chla.shape[0], len(CC.chl_grid_x), len(CC.chl_grid_y)))

    # Interpolate the chla data to this grid for each time step
    for t in range(chla.shape[0]):
        # print("Interpolating time step: ", t + 1, " of ", chla.shape[0])
        interpolated_chla[t] = griddata(
            (x_cropped, y_cropped), chla_cropped[t], (xx, yy),
            method='linear').reshape(len(CC.chl_grid_x), len(CC.chl_grid_y))

    return interpolated_chla


def crop_data(CC, interpolated_chla):
    """ 
    Crop the interpolated data to the extents of the cell 
    :param CC: Chla Object
    :param interpolated_chla: interpolated chla data
    :return: interpolated_chla_cropped
    """
    minx, miny, maxx, maxy = CC.extents
    interpolated_chla_cropped = interpolated_chla[:, (
        CC.chl_grid_y >= miny) & (CC.chl_grid_y <= maxy), :][:, :, (
            CC.chl_grid_x >= minx) & (CC.chl_grid_x <= maxx)]

    CC.chl_grid_x = CC.chl_grid_x[(CC.chl_grid_x >= minx)
                                  & (CC.chl_grid_x <= maxx)]

    CC.chl_grid_y = CC.chl_grid_y[(CC.chl_grid_y >= miny)
                                  & (CC.chl_grid_y <= maxy)]

    interpolated_chla_cropped[interpolated_chla_cropped < 0] = np.nan
    return interpolated_chla_cropped


def store_interpolated_data(CC, file, interpolated_chla, time, file_path):
    """
    Store the interpolated data in a netcdf file
    :param CC: Chla Object
    :param file: file path to the original netcdf file
    :param interpolated_chla: interpolated chla data
    :param time: time
    """
    # Assuming you have your data loaded as interpolated_chla for 3413 coordinates x and y

    # Create an xarray dataset
    ds = xr.Dataset(
        {
            "chla": (["time", "x", "y"], interpolated_chla),
        },
        coords={
            "x": CC.chl_grid_x,
            "y": CC.chl_grid_y,
            # "time": np.arange(interpolated_chla.shape[0]),
            "time": time,
        },
    )

    # Save the dataset to a NetCDF file
    # file_path = os.path.join(CC.data_folder, "interpolated_data",
    #                          CC.short_name, str(CC.active_cell))

    file_name = file.split('/')[-1].split('.')[0]

    # if file_path does not exist, create it
    if not os.path.exists(file_path):
        os.makedirs(file_path)

    file = os.path.join(file_path, file_name + ".nc")

    print("Saving data to file: ", file)
    ds.to_netcdf(file)
    # Print a confirmation message
    print("Data saved to file")


def interpolate_files(CC, source_folder):
    file_list = get_file_list(source_folder)

    lons, lats = get_lat_lon(file_list[0])
    x, y = reproject(lons, lats)

    for file in file_list:
        print("Working on file: ", file)
        file_name = file.split('/')[-1].split('.')[0]

        print("Getting chla stack...")
        chla, time = load_file_xr(file)

        for cell in CC.cells:
            CC.active_cell = cell

            print("Interpolating cell: ", cell)

            cell_path = os.path.join(CC.data_folder, "interpolated_data",
                                     CC.short_name, cell)

            # if cell_path does not exist, create it
            if not os.path.exists(cell_path):
                os.makedirs(cell_path)

            interp_file = os.path.join(cell_path, file_name + ".nc")
            print("Working on :", interp_file)

            # Check if interp_file file exists. If it does, skip this cell
            if os.path.exists(interp_file):
                print(interp_file, "File already exists. Skipping cell...")
                continue

            CC.extents = grid.grid_cell_ID_to_cell_bounds(
                CC.active_cell, CC.grid_dict)
            buffer_extents_and_create_grid(CC)

            print("Interpolating data...")
            interpolated_chla = interpolate_data(CC, x, y, chla)

            print("Cropping data...")
            interpolated_chla = crop_data(CC, interpolated_chla)
            print("Storing data...")
            store_interpolated_data(CC, file, interpolated_chla, time,
                                    cell_path)
            print("Done! Woohoo!")

            del interpolated_chla

        del chla, time

    return


def get_file_list(source_folder):
    # get a list of files in data_folder
    file_list = glob.glob(source_folder + '/*.nc')
    file_list.sort()
    return file_list


def inperpolate_cells(CC, source_folder):
    file_list = get_file_list(source_folder)

    for cell in CC.cells:
        print("Interpolating cell: ", cell)
        CC.active_cell = cell
        interpolate_files(CC, file_list)
    return
