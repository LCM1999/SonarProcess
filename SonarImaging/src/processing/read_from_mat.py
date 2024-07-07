import numpy as np
import vtk, vtkmodules
from vtkmodules.util import numpy_support
import scipy
import cv2
import config

'''
def get_diffuse_grids(file_path):
    data = scipy.io.loadmat(file_path)
    diffuse = data['diffuse'][0][0]
    grids = diffuse['grid'][0]
    return grids

def get_diffuse_data(file_path, value_type):
    grids = get_diffuse_grids(file_path)
    for grid in grids:
        if grid[0][0]['type'][0] == value_type:
            return grid[0][0]['v']
    return np.array([[]])
'''

######################################################################################################

def get_imaging(file_path, value_type):
    data = scipy.io.loadmat(file_path)
    covis = data['imaging'][0][0]
    grid = covis['grid'][0][0]
    config.imaging['data'] = grid[value_type].transpose(1,0,2)
    config.imaging['data_cut'] = config.imaging['data']
    config.imaging['bounds'] = grid['axis'][0]
    config.imaging['cut_bounds'] = grid['axis'][0]
    spacing = grid['spacing'][0][0]
    config.imaging['spacing'] = [spacing['dx'][0][0],spacing['dy'][0][0],spacing['dz'][0][0]]


def get_doppler_data(file_path):
    ...


def get_bathy_from_diffuse(file_path):
    data = scipy.io.loadmat(file_path)
    diffuse = data['diffuse'][0][0]
    grids = diffuse['grid'][0]
    for grid in grids:
        if grid[0][0]['type'][0] == 'bathy_filt':
            config.bathy['data'] = grid[0][0]['v'].transpose()
            config.bathy['bounds'] = grid[0][0]['axis'][0,0:4]
            config.bathy['spacing'] = [grid[0][0]['spacing'][0]['dx'][0][0][0],grid[0][0]['spacing'][0]['dy'][0][0][0]]


#########################################################################################################


'''
def get_bathymetry(file_path):
    data = scipy.io.loadmat(file_path)
    covis = data['covis'][0][0]
    grid = covis['grid'][0][0]
    value = grid['v']
    return value

'''
