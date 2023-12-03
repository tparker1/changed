import numpy as np


def create_grid_from_extents(extents, posting, epsg):

    # the grid is created using the grid frame from bedmachine
    minBedMachineX = -652925
    maxBedMachineX = 879625
    minBedMachineY = -3384425
    maxBedMachineY = -632675

    x = np.arange(minBedMachineX, maxBedMachineX, posting)
    y = np.arange(minBedMachineY, maxBedMachineY, posting)

    x = x[x >= extents[0]]
    x = x[x <= extents[2]]

    y = y[y >= extents[1]]
    y = y[y <= extents[3]]
    return (x, y)


def create_grid_from_extents(IC_object):
    import toolbox.grid_generation as grid

    # get the bounds of the entire ice sheet
    bbox = grid.grid_bounds(
        IC_object.icesheet_name)  # min_x, min_y, max_x, max_y

    # generate the grid
    x = np.arange(bbox[0], bbox[2], IC_object.posting)
    y = np.arange(bbox[1], bbox[3], IC_object.posting)

    # crop the grid to the cell extents
    x = x[x >= IC_object.extents[0]]
    x = x[x < IC_object.extents[2]]

    y = y[y >= IC_object.extents[1]]
    y = y[y < IC_object.extents[3]]

    if IC_object.data_type.lower() == 'velocity':
        IC_object.velocity_grid_x, IC_object.velocity_grid_y = x, y
    elif IC_object.data_type.lower() == 'elevation':
        IC_object.elevation_grid_x, IC_object.elevation_grid_y = x, y
    elif IC_object.data_type.lower() == 'chlorophyll':
        IC_object.chl_grid_x, IC_object.chl_grid_y = x, y

    else:
        print(
            'Warning: Grid not created. \nSet the object data type to either "velocity" or "elevation"'
        )

    return


"""
def create_velocity_grids_from_extents(GD_object):
    velocity_grid_x, velocity_grid_y = \
        create_grid_from_extents(GD_object.extents, GD_object.velocity_grid_posting, GD_object.velocity_grid_epsg)

    GD_object.velocity_grid_x = velocity_grid_x
    GD_object.velocity_grid_y = velocity_grid_y

def create_elevation_grids_from_extents(GD_object):
    elevation_grid_x, elevation_grid_y = \
        create_grid_from_extents(GD_object.extents, GD_object.elevation_grid_posting, GD_object.elevation_grid_epsg)

    GD_object.elevation_grid_y = elevation_grid_y
    GD_object.elevation_grid_x = elevation_grid_x
"""
