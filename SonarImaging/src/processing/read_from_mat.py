import numpy as np
import vtk, vtkmodules
from vtkmodules.util import numpy_support
import scipy
import cv2
# import SonarImaging.src.config as config

from SonarProcess.Imaging.buildImageGrid import *

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

# def get_imaging(file_path, value_type):
#     data = scipy.io.loadmat(file_path)
#     covis = data['imaging'][0][0]
#     grid = covis['grid'][0][0]
#     config.imaging['data'] = grid[value_type].transpose(1, 0, 2)
#     config.imaging['data_cut'] = config.imaging['data']
#     config.imaging['bounds'] = grid['axis'][0]
#     config.imaging['cut_bounds'] = grid['axis'][0]
#     spacing = grid['spacing'][0][0]
#     config.imaging['spacing'] = [spacing['dx'][0][0], spacing['dy'][0][0], spacing['dz'][0][0]]

def get_imaging_image(file_path):
    data = scipy.io.loadmat(file_path)
    imaging = {}
    data = data['covis']['grid'][0][0][0][0]
    imaging['units'] = [data[2][0][0][i][0] for i in range(2)]
    imaging['bounds'] = [data[4][0][0][i][0][0] for i in range(6)]
    imaging['spacing'] = [data[5][0][0][i][0][0] for i in range(3)]
    imaging['x'] = data[6]
    imaging['y'] = data[7]
    imaging['z'] = data[8]
    imaging['v'] = data[9]
    imaging['w'] = data[10]
    imaging['shape'] = [data[11][0][i] for i in range(3)]
    imaging['name'] = data[13][0]

    return imaging

def get_doppler_data(file_path):
    data = scipy.io.loadmat(file_path)
    doppler = {}
    data = data['covis']['grid'][0][0][0][0][0][0]

    doppler['unit'] = [data[2][0][0][i][0] for i in range(2)]
    doppler['bounds'] = [data[4][0][0][i][0][0] for i in range(6)]
    doppler['spacing'] = [data[5][0][0][i][0][0] for i in range(3)]
    doppler['x'] = np.transpose(data[6], (1, 0, 2))
    doppler['y'] = np.transpose(data[7], (1, 0, 2))
    doppler['z'] = np.transpose(data[8], (1, 0, 2))
    doppler['v'] = np.transpose(data[9], (1, 0, 2))
    doppler['v_filt'] = np.transpose(data[10], (1, 0, 2))
    doppler['w'] = np.transpose(data[11], (1, 0, 2))
    doppler['std'] = np.transpose(data[12], (1, 0, 2))
    doppler['vr'] = np.transpose(data[13], (1, 0, 2))
    doppler['covar'] = np.transpose(data[14], (1, 0, 2))
    doppler['shape'] = [data[15][0][1], data[15][0][0], data[15][0][2]]
    doppler['name'] = data[17][0]
    # doppler['offset_covar'] = [data[18][0][i] for i in range(2)]

    return doppler


def get_doppler_image(file_path):
    data = scipy.io.loadmat(file_path)
    doppler = {}
    data = data['covis'][0][0][5][0][1][0][0]

    doppler['unit'] = [data[2][0][0][i][0] for i in range(2)]
    doppler['bounds'] = [data[4][0][0][i][0][0] for i in range(6)]
    doppler['spacing'] = [data[5][0][0][i][0][0] for i in range(3)]
    doppler['x'] = np.transpose(data[6], (1, 0, 2))
    doppler['y'] = np.transpose(data[7], (1, 0, 2))
    doppler['z'] = np.transpose(data[8], (1, 0, 2))
    doppler['v'] = np.transpose(data[9], (1, 0, 2))
    doppler['v_filt'] = np.transpose(data[10], (1, 0, 2))
    doppler['w'] = np.transpose(data[11], (1, 0, 2))
    doppler['std'] = np.transpose(data[12], (1, 0, 2))
    doppler['vr'] = np.transpose(data[13], (1, 0, 2))
    doppler['covar'] = np.transpose(data[14], (1, 0, 2))
    doppler['shape'] = [data[15][0][1], data[15][0][0], data[15][0][2]]
    doppler['name'] = data[17][0]

    return doppler


# def get_bathy_from_diffuse(file_path):
#     data = scipy.io.loadmat(file_path)
#     diffuse = data['diffuse'][0][0]
#     grids = diffuse['grid'][0]
#     for grid in grids:
#         if grid[0][0]['type'][0] == 'bathy_filt':
#             config.bathy['data'] = grid[0][0]['v'].transpose()
#             config.bathy['bounds'] = grid[0][0]['axis'][0, 0:4]
#             config.bathy['spacing'] = [grid[0][0]['spacing'][0]['dx'][0][0][0], grid[0][0]['spacing'][0]['dy'][0][0][0]]


#########################################################################################################


'''
def get_bathymetry(file_path):
    data = scipy.io.loadmat(file_path)
    covis = data['covis'][0][0]
    grid = covis['grid'][0][0]
    value = grid['v']
    return value

'''

if __name__ == "__main__":
    doppler_dict = get_doppler_data("F:/Covis/pythonProject/Inputs/APLUWCOVISMBSONAR001_20141006T210836.397Z-DOPPLER.mat")
    doppler_image = build_image_grid(array_dict={
        'v': doppler_dict['v'],
        'v_filt': doppler_dict['v_filt'],
        'w': doppler_dict['w'],
        'std': doppler_dict['std'],
        'vr': doppler_dict['vr'],
        'covar': doppler_dict['covar']
    }, bounds=doppler_dict['bounds'], spacing=doppler_dict['spacing'])
    save_image(doppler_image, "F:/data/", 'doppler.vti')

    intensity = get_doppler_image("F:/Covis/pythonProject/Inputs/APLUWCOVISMBSONAR001_20141006T210836.397Z-DOPPLER.mat")
    intensity_image = build_image_grid(array_dict={
        'v': intensity['v'],
        'v_filt': intensity['v_filt'],
        'w': intensity['w'],
        'std': intensity['std'],
        'vr': intensity['vr'],
        'covar': intensity['covar']
    }, bounds=intensity['bounds'], spacing=intensity['spacing'])

    save_image(intensity_image, "F:/data/", 'intensity.vti')
    # imaging = get_imaging_image("F:/Covis/pythonProject/Inputs/APLUWCOVISMBSONAR001_20141006T210021.023Z-IMAGING.mat")
    # image = build_image_grid(array_dict={
    #     'v': imaging['v'],
    #     'w': imaging['w'],
    # }, bounds=imaging['bounds'], spacing=imaging['spacing'])
    #
    # save_image(image, "F:/data/", "imaging.vti")
