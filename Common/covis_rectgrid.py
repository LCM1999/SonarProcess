import numpy as np

def covis_rectgrid(grd):
    eps = 0.00001

    xmin = grd['bounds']['xmin'] if 'xmin' in grd['bounds'] else 0
    xmax = grd['bounds']['xmax'] if 'xmax' in grd['bounds'] else 0
    ymin = grd['bounds']['ymin'] if 'ymin' in grd['bounds'] else 0
    ymax = grd['bounds']['ymax'] if 'ymax' in grd['bounds'] else 0
    zmin = grd['bounds']['zmin'] if 'zmin' in grd['bounds'] else 0
    zmax = grd['bounds']['zmax'] if 'zmax' in grd['bounds'] else 0

    dx = grd['spacing']['dx'] if 'dx' in grd['spacing'] else 0
    dy = grd['spacing']['dy'] if 'dy' in grd['spacing'] else 0
    dz = grd['spacing']['dz'] if 'dz' in grd['spacing'] else 0

    dims = grd['dimensions'] if 'dimensions' in grd else 3

    x = np.arange(xmin, xmax + eps * dx, dx) if dx != 0 else np.array([xmin])
    y = np.arange(ymin, ymax + eps * dy, dy) if dy != 0 else np.array([ymin])
    z = np.arange(zmin, zmax + eps * dz, dz) if dz != 0 else np.array([zmin])

    # use meshgrid to define the mesh
    if dims == 3:
        grd['x'], grd['y'], grd['z'] =np. meshgrid(x, y, z)    # define grid mesh coords
    if dims == 2:
        grd['x'], grd['y'] = np.meshgrid(x,y)    # define grid mesh coords

    # initialize the grid value and weight matrix
    grd['v'] = np.zeros(grd['x'].shape)    # define an empty grid value matrix
    grd['w'] = np.zeros(grd['x'].shape)    # define an empty grid weight matrix

    # set the bounds and size
    grd['size'] = grd['v'].shape
    grd['axis'] = [xmin, xmax, ymin, ymax, zmin, zmax]

    return grd
