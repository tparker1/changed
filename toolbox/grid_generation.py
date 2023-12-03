import shapefile
from osgeo import osr
import os
import xarray as xr
import numpy as np


# Function to generate a grid bound dictionary for a given region (greenland or antarctic)
def grid_bounds(region):
    # min_x, min_y, max_x, max_y
    grid_dict = {
        "Greenland": [-652925, -3384425, 850000, -594400],
        "Antarctic": [-3333000, -3333000, 3333000, 3333000]
    }
    return grid_dict[region.title()]


def generate_grid_dictionary(min_x, min_y, max_x, max_y, resolution_km):

    step = resolution_km * 1000  # Convert to meters

    grid_dict = {}
    for row in range(0, int((max_y - min_y) / step) + 1):
        y_1 = min_y + (row * step)
        y_2 = min_y + ((row + 1) * step)
        for col in range(0, int((max_x - min_x) / step) + 1):
            x_1 = min_x + (col * step)
            x_2 = min_x + ((col + 1) * step)
            grid_dict['r_' + str(row).zfill(2) + '_c_' +
                      str(col).zfill(2)] = [x_1, y_1, x_2, y_2]

    return grid_dict


def create_shapefile(IC):

    if IC.data_type.lower() == 'velocity':
        epsg = IC.velocity_grid_epsg
        resolution = IC.velocity_grid_posting
    elif IC.data_type.lower() == 'elevation':
        epsg = IC.elevation_grid_epsg
        resolution = IC.elevation_grid_posting
    elif IC.data_type.lower() == 'chlorophyll':
        epsg = IC.chl_grid_epsg
        resolution = IC.chl_grid_posting

    else:
        print(
            'Warning: Grid not created. \nSet the object data type to either "Velocity" or "Elevation"'
        )
        return

    osr.UseExceptions()
    output_path = os.path.join(IC.project_folder, "Metadata", IC.icesheet_name,
                               IC.data_type.title(),
                               "Grid_" + str(resolution) + "km")

    fields = ['cell_id']
    fieldTypes = ['C']
    records = []
    polygons = []

    col_range = int(list(IC.grid_dict.keys())[-1].split('_')[-1]) + 1
    row_range = int(list(IC.grid_dict.keys())[-1].split('_')[1]) + 1

    for col in range(col_range):
        for row in range(row_range):
            cell_ID = 'r_' + '{:02d}'.format(row) + '_c_' + '{:02d}'.format(
                col)
            bounds = IC.grid_dict[cell_ID]
            polygon = [[bounds[0], bounds[1]], [bounds[2], bounds[1]],
                       [bounds[2], bounds[3]], [bounds[0], bounds[3]]]
            records.append(cell_ID)
            polygons.append(polygon)

    w = shapefile.Writer(output_path)

    for ff in range(len(fields)):
        w.field(fields[ff], fieldTypes[ff])

    for p in range(len(polygons)):
        w.record(records[p])
        w.poly([polygons[p]])

    w.close()

    spatialRef = osr.SpatialReference()
    spatialRef.ImportFromEPSG(IC.epsg)
    spatialRef.MorphToESRI()
    file = open(output_path + '.prj', 'w')
    file.write(spatialRef.ExportToWkt())
    file.close()
    print("Shapefile saved to " + output_path)
    return


def grid_cell_ID_to_cell_bounds(cell, grid_dict):
    return grid_dict[cell]


def create_grid(IC):
    bbox = grid_bounds(IC.icesheet_name.title())

    if IC.data_type.lower() == 'velocity':
        IC.grid_dict = generate_grid_dictionary(bbox[0], bbox[1], bbox[2],
                                                bbox[3],
                                                IC.velocity_grid_posting)
        print("Velocity grid dictionary generated for " +
              str(IC.velocity_grid_posting) + "km resolution")
        create_shapefile(IC)

    elif IC.data_type.lower() == 'elevation':
        IC.grid_dict = generate_grid_dictionary(bbox[0], bbox[1], bbox[2],
                                                bbox[3],
                                                IC.elevation_grid_posting)
        print("Elevation grid dictionary generated for " +
              str(IC.elevation_grid_posting) + "km resolution")
        create_shapefile(IC)

    elif IC.data_type.lower() == 'chlorophyll':
        IC.grid_dict = generate_grid_dictionary(bbox[0], bbox[1], bbox[2],
                                                bbox[3], IC.chl_grid_posting)
        print("Chlorophyll grid dictionary generated for " +
              str(IC.chl_grid_posting) + "km resolution")
        create_shapefile(IC)

    else:
        print(
            'Warning: Grid not created. \ndata_type "', IC.data_type.lower(),
            '" not supported. \nSet the object data type to either "Velocity", "Elevation", or "Chlorophyll"'
        )
    return


# Store the cells of interest as IC object attribute (IC.cells)
def store_cells(IC, grid):
    IC.cells = []
    for row in range(grid[0], grid[2] + 1):
        for col in range(grid[1], grid[3] + 1):
            cell_ID = 'r_' + '{:02d}'.format(row) + '_c_' + '{:02d}'.format(
                col)
            if cell_ID not in IC.cells:
                IC.cells.append(cell_ID)
    IC.cells = sorted(IC.cells)
    print("Cells: ", IC.cells)
    return


def organize_quilt_grid_input(cell_IDs):
    lowest_row = 1e10
    highest_row = -1e10
    lowest_col = 1e10
    highest_col = -1e10
    for cell_ID in cell_IDs:
        row = int(cell_ID.split('_')[1])
        col = int(cell_ID.split('_')[3])
        if row < lowest_row:
            lowest_row = row
        if row > highest_row:
            highest_row = row
        if col < lowest_col:
            lowest_col = col
        if col > highest_col:
            highest_col = col
    row_col_extents = [lowest_row, highest_row, lowest_col, highest_col]

    ordered_cell_ID_grid = []
    for row in range(lowest_row, highest_row + 1):
        row_cells = []
        for col in range(lowest_col, highest_col + 1):
            cell_test = 'r_' + '{:02d}'.format(row) + '_c_' + '{:02d}'.format(
                col)
            if cell_test in cell_IDs:
                row_cells.append(cell_test)
            else:
                row_cells.append('')
        ordered_cell_ID_grid.append(row_cells)
    return (ordered_cell_ID_grid, row_col_extents)


def create_input_file_path_dict(IC):
    output_path = os.path.join(
        IC.data_folder, 'Output',
        IC.icesheet_name)  # + IC.glacier_name + "_stack.nc")
    output_files = []
    input_file_path_dict = {}

    for cell in IC.cells:
        output_file = os.path.join(output_path, cell,
                                   IC.short_name + "_stack.nc")
        output_files.append(output_file)

        #check if file exists
        if os.path.exists(output_file):
            input_file_path_dict[cell] = output_file
        else:
            print("WARNING: file does not exist: ", output_file)
    #print(input_file_path_dict)

    return input_file_path_dict, os.path.join(output_path,
                                              IC.data_type.title())


def quilt_grids_by_path(IC_object,
                        input_file_paths,
                        grid_variable_label,
                        resolution,
                        no_data_value=0):
    lowest_row = 1e10
    highest_row = -1e10
    lowest_col = 1e10
    highest_col = 1e10
    cell_IDs, grid_dict = IC_object.cells, IC_object.grid_dict

    for cell_ID in cell_IDs:
        row = int(cell_ID.split('_')[1])
        col = int(cell_ID.split('_')[3])
        if row < lowest_row:
            lowest_row = row
        if row >= highest_row:  # >= to make sure the highest row/col is included if there is only 1 row/col
            highest_row = row
        if col < lowest_col:
            lowest_col = col
        if col >= lowest_col:
            highest_col = col

    print("Lowest row: ", lowest_row)
    print("Highest row: ", highest_row)
    print("Lowest col: ", lowest_col)
    print("Highest col: ", highest_col)

    UR_cell = 'r_' + '{:02d}'.format(highest_row) + '_c_' + '{:02d}'.format(
        highest_col)
    UR_bounds = grid_cell_ID_to_cell_bounds(UR_cell, grid_dict)
    LL_cell = 'r_' + '{:02d}'.format(lowest_row) + '_c_' + '{:02d}'.format(
        lowest_col)
    LL_bounds = grid_cell_ID_to_cell_bounds(LL_cell, grid_dict)

    print("Resolution in function: ", resolution)
    x = np.arange(LL_bounds[0], UR_bounds[2], resolution)
    y = np.arange(LL_bounds[1], UR_bounds[3], resolution)
    print("LL bounds: ", LL_bounds)
    print("UR bounds: ", UR_bounds)

    file_path = input_file_paths[cell_IDs[0]]
    ds = xr.open_dataset(file_path)
    grid = np.array(ds.data_vars[grid_variable_label])

    print("Grid shape: ", np.shape(grid))
    print("Grid[0] shape: ", np.shape(grid[0]))
    time_len = np.shape(grid)[0]

    if IC_object.short_name == 'MEaSUREs Greenland Monthly Velocity':
        t = ds.variables['time_start'][:]
    else:
        t = ds.variables['time'][:]

    if no_data_value != 'nan':
        main_grid = no_data_value * np.ones((time_len, len(y), len(x)))
    else:
        main_grid = np.zeros((time_len, len(y), len(x)))
        main_grid[::] = np.nan

    for cell in input_file_paths:
        print("Loading cell ", cell)
        file_path = input_file_paths[cell]
        cell_ID = cell
        print("File path: ", file_path)

        ds = xr.open_dataset(file_path)
        grid = np.array(ds.data_vars[grid_variable_label])
        ds = None

        row = int(cell_ID.split('_')[1])
        col = int(cell_ID.split('_')[3])

        jG = (row - lowest_row) * np.shape(grid)[2]
        iG = (col - lowest_col) * np.shape(grid)[1]

        main_grid[:, jG:jG + np.shape(grid)[1],
                  iG:iG + np.shape(grid)[2]] = grid

    #main_grid[np.isnan(main_grid)] = 0

    return (t, x, y, main_grid)


def quilt_grids_and_output(IC, grid_variable_label, no_data_value=0):
    import toolbox.store_nc as store

    input_file_paths, output_path = create_input_file_path_dict(IC)
    file_name = IC.nickname + "_stack.nc"

    if os.path.exists(os.path.join(output_path, file_name)):
        print("WARNING: file already exists: ",
              os.path.join(output_path, file_name))
        print(
            "Give the object a unique nickname or delete the existing file and try again."
        )
        return

    output_file = os.path.join(output_path, file_name)
    if IC.data_type.lower() == 'velocity':
        resolution = IC.velocity_grid_posting
    elif IC.data_type.lower() == 'elevation':
        resolution = IC.elevation_grid_posting
    else:
        print(
            'Warning: Grid not created. \nSet the object data type to either "velocity" or "elevation"'
        )
        return

    print("Resolution in quilt_grids_and_output: ", resolution)
    #if IC.short_name == 'ATL15 Antarctic Elevation':
    t, x, y, main_grid = quilt_grids_by_path(IC, input_file_paths,
                                             grid_variable_label, resolution,
                                             no_data_value)
    store.output_grid_with_geo_reference_timedim_atl(IC, output_file, t, x, y,
                                                     main_grid,
                                                     grid_variable_label,
                                                     resolution)
    #elif IC.short_name == 'MEaSUREs Greenland Monthly Velocity':
    #    t, x, y, main_grid = quilt_grids_by_path(IC, input_file_paths, grid_variable_label, resolution, no_data_value)
    #    store.output_grid_with_geo_reference_timedim_measures(IC, output_file, t, x, y, main_grid, grid_variable_label, resolution)

    return (main_grid)
