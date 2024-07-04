import vtk, vtkmodules
import numpy as np
import joblib
import cv2
from vtkmodules.util import numpy_support, vtkImageImportFromArray
import scipy
from processing import read_from_mat, feature_detection
import time
import config


def to_vtk_data(data):
    vtk_data = vtk.vtkImageData()
    vtk_data.SetDimensions(data.shape) 
    vtk_data.SetSpacing(config.imaging['spacing'][0], config.imaging['spacing'][1], config.imaging['spacing'][2])
    vtk_data.SetOrigin(config.imaging['cut_bounds'][0], config.imaging['cut_bounds'][2], config.imaging['cut_bounds'][4])
 
    data = data.astype(np.float32)
    data = data.transpose(2,1,0)
    vtk_array = numpy_support.numpy_to_vtk(data.ravel(), deep=True, array_type=vtk.VTK_FLOAT)
    vtk_array.SetNumberOfComponents(1)
    vtk_array.SetName("data")

    vtk_data.GetPointData().SetScalars(vtk_array)
    return vtk_data



def imaging_cutting(bounds=[-30,20,-30,20,0,20]):
    config.imaging['cut_bounds'] = bounds
    dx1 = int((bounds[0]-config.imaging['bounds'][0])/config.imaging['spacing'][0])
    dx2 = int((bounds[1]-config.imaging['bounds'][1])/config.imaging['spacing'][0])
    dy1 = int((bounds[2]-config.imaging['bounds'][2])/config.imaging['spacing'][1])
    dy2 = int((bounds[3]-config.imaging['bounds'][3])/config.imaging['spacing'][1])
    sz = config.imaging['data'].shape
    config.imaging['data_cut'] = config.imaging['data'][dx1:sz[0]-dx2,dy1:sz[1]-dy2,:]
    config.imaging['vtk_data'] = to_vtk_data(config.imaging['data_cut'])

